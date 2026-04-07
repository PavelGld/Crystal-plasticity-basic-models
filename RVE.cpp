#include "RVE.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <unsupported/Eigen/MatrixFunctions>

#ifdef _OPENMP
#include <omp.h>
#endif

RVE::RVE(const std::string& config_file) {
    std::ifstream file(config_file);

    int num_materials;
    file >> num_materials;

    materials.resize(num_materials);
    for (int m = 0; m < num_materials; m++) {
        std::string mat_file;
        int num_grains_of_material;
        file >> mat_file >> num_grains_of_material;

        auto mat = std::make_shared<MaterialParams>(mat_file);
        materials[m] = mat;

        for (int g = 0; g < num_grains_of_material; g++) {
            Eigen::Quaterniond random_q = Eigen::Quaterniond::UnitRandom();
            grains.emplace_back(mat, random_q);
        }
    }
    
    stress_av.setZero();
    W_av.setZero();
    Dp_av.setZero();
    F_av.setIdentity();
    A_affine.setIdentity();
    time = 0;
    acc_strain = 0;
    
    total_volume = 0;
    for (const auto& grain : grains) {
        total_volume += grain.get_volume();
    }
}

void RVE::Simulate_Step(double dt, const Eigen::Matrix3d& L) {
    stress_av.setZero();
    W_av.setZero();
    Dp_av.setZero();
    Dp_v_buf.setZero();
    
    int n = static_cast<int>(grains.size());
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n; i++) {
        grains[i].Simulate_Step(dt, L);
    }
    
    for (int i = 0; i < n; i++) {
        double weight = grains[i].get_volume() / total_volume;
        temp_3x3_buf.noalias() = grains[i].get_stress_lab();
        stress_av.noalias() += temp_3x3_buf * weight;

        // Осредняем Dp для диагностики
        Dp_v_buf += weight * grains[i].get_Dp_v();
    }

    Dp_av(0,0) = Dp_v_buf(0);
    Dp_av(1,1) = Dp_v_buf(1);
    Dp_av(2,2) = Dp_v_buf(2);
    Dp_av(1,2) = Dp_av(2,1) = Dp_v_buf(3) / 2;
    Dp_av(0,2) = Dp_av(2,0) = Dp_v_buf(4) / 2;
    Dp_av(0,1) = Dp_av(1,0) = Dp_v_buf(5) / 2;

    // acc_strain из полной деформации (для контроля цикла)
    D_buf = 0.5 * (L + L.transpose());
    double tr_D = D_buf.trace();
    D_buf.diagonal().array() -= tr_D / 3.0;
    double D_eq = std::sqrt(2.0 / 3.0 * D_buf.squaredNorm());
    acc_strain += D_eq * dt;

    time += dt;
    
    exp_buf.noalias() = dt * L;
    F_av.noalias() = exp_buf.exp() * F_av;
}

void RVE::Simulate_Step_Soft(double dt, const Eigen::Matrix3d& dS_dt) {
    int n = static_cast<int>(grains.size());
    
    // 1. Вычислить среднюю упругость (в ЛСК, Фойгт)
    P_av_buf.setZero();
    
    for (int i = 0; i < n; i++) {
        double w = grains[i].get_volume() / total_volume;
        grains[i].get_elasticity_lab_into(temp_P_buf);
        P_av_buf.noalias() += w * temp_P_buf;
    }
    
    // 2. Вычислить D из формулы: D = P^-1 : dS/dt + Dp
    dS_dt_v_buf << dS_dt(0,0), dS_dt(1,1), dS_dt(2,2), 
                   dS_dt(1,2)+dS_dt(2,1), dS_dt(0,2)+dS_dt(2,0), 
                   dS_dt(0,1)+dS_dt(1,0);
    
    // LU-инверсия без аллокаций
    lu_buf.compute(P_av_buf);
    D_v_buf.noalias() = lu_buf.solve(dS_dt_v_buf);
    D_v_buf.noalias() += Dp_av_old_v;
    
    // 3. Преобразовать в L (спин с прошлого шага)
    D_buf << D_v_buf(0), D_v_buf(5)/2, D_v_buf(4)/2,
             D_v_buf(5)/2, D_v_buf(1), D_v_buf(3)/2,
             D_v_buf(4)/2, D_v_buf(3)/2, D_v_buf(2);
    
    L_buf.noalias() = D_buf + W_av_old;
    
    // 4. Параллельный расчёт зёрен
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n; i++) {
        grains[i].Simulate_Step(dt, L_buf);
    }
    
    // 5. Осреднить результаты
    stress_v_buf.setZero();
    Dp_v_buf.setZero();
    W_sum_buf.setZero();
    
    for (int i = 0; i < n; i++) {
        double w = grains[i].get_volume() / total_volume;
        
        // tensorToVoigt inline для get_stress_lab()
        const Eigen::Matrix3d& s = grains[i].get_stress_lab();
        stress_v_buf(0) += w * s(0,0);
        stress_v_buf(1) += w * s(1,1);
        stress_v_buf(2) += w * s(2,2);
        stress_v_buf(3) += w * (s(1,2) + s(2,1));
        stress_v_buf(4) += w * (s(0,2) + s(2,0));
        stress_v_buf(5) += w * (s(0,1) + s(1,0));
        
        Dp_v_buf.noalias() += w * grains[i].get_Dp_v();
        
        // W_sum += w * 0.5 * (L - L^T)
        W_sum_buf(0,1) += w * 0.5 * (L_buf(0,1) - L_buf(1,0));
        W_sum_buf(0,2) += w * 0.5 * (L_buf(0,2) - L_buf(2,0));
        W_sum_buf(1,2) += w * 0.5 * (L_buf(1,2) - L_buf(2,1));
        W_sum_buf(1,0) = -W_sum_buf(0,1);
        W_sum_buf(2,0) = -W_sum_buf(0,2);
        W_sum_buf(2,1) = -W_sum_buf(1,2);
    }
    
    // 6. Сохранить результаты
    stress_av(0,0) = stress_v_buf(0);
    stress_av(1,1) = stress_v_buf(1);
    stress_av(2,2) = stress_v_buf(2);
    stress_av(1,2) = stress_av(2,1) = stress_v_buf(3) / 2;
    stress_av(0,2) = stress_av(2,0) = stress_v_buf(4) / 2;
    stress_av(0,1) = stress_av(1,0) = stress_v_buf(5) / 2;
    
    stress_av_old_v = stress_v_buf;
    Dp_av_old_v = Dp_v_buf;
    W_av_old = W_sum_buf;
    W_av = W_sum_buf;
    
    // Dp_av
    Dp_av(0,0) = Dp_v_buf(0);
    Dp_av(1,1) = Dp_v_buf(1);
    Dp_av(2,2) = Dp_v_buf(2);
    Dp_av(1,2) = Dp_av(2,1) = Dp_v_buf(3) / 2;
    Dp_av(0,2) = Dp_av(2,0) = Dp_v_buf(4) / 2;
    Dp_av(0,1) = Dp_av(1,0) = Dp_v_buf(5) / 2;
    
    time += dt;
    
    exp_buf.noalias() = dt * L_buf;
    F_av.noalias() = exp_buf.exp() * F_av;
}
