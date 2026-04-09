Баги
1. get_elasticity_lab() уничтожает свой кэш (Grain.cpp)

const Grain::Matrix6d& Grain::get_elasticity_lab() const {
    if (!T_voigt_valid) {
        computeVoigtRotationInPlace(T_voigt_buf);  // T_voigt_buf = T (матрица поворота)
        T_voigt_valid = true;
    }
    Buf_P.noalias() = T_voigt_buf * params->P_v;
    T_voigt_buf.noalias() = Buf_P * T_voigt_buf.transpose(); // ← ПЕРЕЗАПИСЫВАЕТ T результатом!
    return T_voigt_buf;
}

После первого вызова T_voigt_buf уже не матрица поворота, а T·P·T^T. Второй вызов (если флаг не сброшен) вычислит (T·P·T^T)·P·(T·P·T^T)^T — полную чушь. В отличие от корректной get_elasticity_lab_into, которая пишет в отдельный result. Исправление:

const Grain::Matrix6d& Grain::get_elasticity_lab() const {
    if (!T_voigt_valid) {
        computeVoigtRotationInPlace(T_voigt_buf);
        T_voigt_valid = true;
    }
    // Результат пишем в отдельный буфер, не трогая T_voigt_buf
    Buf_P.noalias() = T_voigt_buf * params->P_v;
    P_lab_buf.noalias() = Buf_P * T_voigt_buf.transpose(); // нужен доп. буфер Matrix6d
    return P_lab_buf;
}

2. computeVoigtRotationInPlace — неверные коэффициенты (Grain.cpp)
Используемое Фойгт-соглашение: для касательных компонент v(3..5) = 2·σ_{ij} (инженерные сдвиги). При таком соглашении матрица Bond T должна удовлетворять:
Тип строки P Тип столбца I Правильный коэфф.В коде
1. нормаль (p==q)нормаль (i==j)R²✅
2. сдвиг (p≠q)сдвиг (i≠j)R_{ip}R_{jq}+R_{iq}R_{jp}✅ (суммируется по (i,j)+(j,i))
3. нормаль (p==q)сдвиг (i≠j)R_{ip}R_{jp}❌ даёт 4·R_{ip}R_{jp}
4. сдвиг (p≠q)нормаль (i==j)2·R_{ip}R_{iq}❌ даёт 0.5·R_{ip}R_{iq}


Оба смешанных случая ошибаются в 4 раза (в разных направлениях). Исправление:

void Grain::computeVoigtRotationInPlace(Matrix6d& T) const {
    T.setZero();
    for (int p = 0; p < 3; ++p)
    for (int q = 0; q < 3; ++q) {
        int P = tensorToVoigtIndex(p, q);
        bool P_shear = (p != q);
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            int I = tensorToVoigtIndex(i, j);
            bool I_shear = (i != j);
            double v = R_buf(i, p) * R_buf(j, q);
            if (!P_shear && !I_shear)       T(I, P) += v;          // норм→норм
            else if (P_shear && I_shear)    T(I, P) += v;          // сдв→сдв (суммируется дважды)
            else if (!P_shear && I_shear)   T(I, P) += v * 0.5;    // норм→сдв
            else                            T(I, P) += v * 2.0;    // сдв→норм
        }
    }
}

⚠️ Этот баг влияет только на Simulate_Step_Soft. Если используешь Simulate_Step (с заданным L) — эффекта нет.


3. Двойная интеграция при переходе упругость→пластика (Grain.cpp)
В check_elastic, когда обнаруживается переход (не was_in_plastic):

// Внутри check_elastic:
Sum_time += buft;       // Sum_time = S + buft
integrate(dT);          // dT = dt-buft — интегрирует ОСТАВШЕЕСЯ время
return false;


Возврат в Simulate_Step:

else {
    find_SS(dT);        // dT — ПОЛНЫЙ шаг, не остаток!
    Sum_time += dT;     // Sum_time = S + buft + (dt-S) = dt + buft > dt
}

find_SS вызывается с полным dT вместо остатка dt - buft, и integrate внутри него снова обновит упрочнение и ориентацию — двойной счёт. Фикс:

// В Grain.h добавить:
double last_elastic_buft = 0.0;

// check_elastic: сохранять buft вместо модификации Sum_time
last_elastic_buft = buft;
// НЕ трогать Sum_time здесь

// Simulate_Step:
if (check_elastic(dT)) {
    Sum_time += dT;
} else {
    double dt_plastic = dT - last_elastic_buft;
    Sum_time += dT;  // полный шаг всё равно засчитываем
    if (dt_plastic > 1e-15)
        find_SS(dt_plastic);
}

Оптимизации
4. Передача Lnew по значению (Grain.h / Grain.cpp)

// БЫЛО:
void Simulate_Step(double dt, Eigen::Matrix3d Lnew);

// СТАЛО:
void Simulate_Step(double dt, const Eigen::Matrix3d& Lnew);

Копирование 3×3 матрицы при каждом вызове — мелочь для 1 зерна, но при 1000 зёрен и тысячах шагов складывается.

5. W_sum_buf бессмысленно вычисляется в цикле по зёрнам (RVE.cpp)
В Simulate_Step_Soft спин зависит только от L_buf, одинакового для всех зёрен:


// БЫЛО — внутри цикла по i, накапливается с весами (но L_buf одинаков):
W_sum_buf(0,1) += w * 0.5 * (L_buf(0,1) - L_buf(1,0));
// ...и т.д.

// СТАЛО — один раз после цикла:
W_sum_buf.setZero();
W_sum_buf(0,1) = 0.5 * (L_buf(0,1) - L_buf(1,0));
W_sum_buf(0,2) = 0.5 * (L_buf(0,2) - L_buf(2,0));
W_sum_buf(1,2) = 0.5 * (L_buf(1,2) - L_buf(2,1));
W_sum_buf(1,0) = -W_sum_buf(0,1);
W_sum_buf(2,0) = -W_sum_buf(0,2);
W_sum_buf(2,1) = -W_sum_buf(1,2);
