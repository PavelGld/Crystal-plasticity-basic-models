#pragma once

#include "trajectory_generator.hpp"  // generateRW(), generateSM(), computeInterval()

// ============================================================
//  Как определить тип и seed траектории по номеру эксперимента
//
//  Соотношение RW:SM = 2:1 (как в статье: 20000 RW + 10000 SM)
//  Первые 2/3 от N  → RW,  оставшиеся 1/3 → SM
//
//  Если хотите другое соотношение — измените rw_count ниже.
// ============================================================

// Заменяет оригинальный while-цикл.
//
// Вызов в вашем коде:
//
//   // было:
//   while (s_sum < s_max) {
//       update_L();
//       for (unsigned k = 0; k < n; k++) model.evolve();
//       output_semax();
//   }
//
//   // стало:
//   run_loading_experiment(model, i, N, dt);
//
// Параметры:
//   model  — ваш объект модели (методы: set_L, evolve, output_semax)
//   i      — номер текущего эксперимента (0-based)
//   N      — общее число экспериментов
//   dt     — шаг по деформации (единицы деформации, напр. 1e-4)

template<typename Model>
void run_loading_experiment(Model&       model,
                            unsigned int i,
                            unsigned int N,
                            double       dt = 1e-4)
{
    // ----------------------------------------------------------
    // 1. Определяем тип траектории и seed
    // ----------------------------------------------------------
    unsigned int rw_count = (N * 2) / 3;  // первые 2/3 — RW, остальные — SM
    bool         is_rw    = (i < rw_count);
    unsigned int seed     = i;             // уникальный seed для каждого эксперимента

    // ----------------------------------------------------------
    // 2. Генерируем последовательность тензоров деформации ε
    //    RW: 1001 точка (включая нулевую)
    //    SM: 1000 точек
    // ----------------------------------------------------------
    std::vector<Matrix3d> eps_seq = is_rw
        ? generateRW(seed, /*n_final=*/1001)
        : generateSM(seed, /*n_final=*/1000);

    // ----------------------------------------------------------
    // 3. Переписанный while-цикл
    //
    //    Оригинал:
    //      while (s_sum < s_max) {
    //          update_L();
    //          for (unsigned k = 0; k < n; k++) model.evolve();
    //          output_semax();
    //      }
    //
    //    Теперь вместо s_sum < s_max — итерация по интервалам траектории.
    //    Вместо update_L() — вычисление L_norm из соседних ε и вызов set_L.
    //    n шагов определяется через длину интервала и dt.
    // ----------------------------------------------------------
    for (size_t seg = 0; seg + 1 < eps_seq.size(); ++seg)
    {
        // --- Аналог: update_L() ---
        //
        // Вычисляем нормированный градиент скорости L и число шагов n
        // для перехода от eps_seq[seg] к eps_seq[seg+1].
        //
        // L_norm * dt  →  за n шагов модель пройдёт ~s_int деформации.
        // Остаток (s_int mod dt) < dt — недохождение меньше одного шага.
        IntervalData iv = computeInterval(eps_seq[seg], eps_seq[seg + 1], dt);

        if (iv.n_steps == 0)
            continue;  // интервал короче одного шага dt — пропускаем

        model.set_L(iv.L_norm);  // <-- вместо update_L()

        // --- Аналог: for (k < n) model.evolve() ---
        for (int k = 0; k < iv.n_steps; ++k)
            model.evolve();

        // --- Аналог: output_semax() ---
        model.output_semax();
    }
}
