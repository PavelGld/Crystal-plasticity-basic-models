#include "RVE.h"
#include <iostream>
#include <fstream>
#include <unsupported/Eigen/MatrixFunctions>

#ifdef _OPENMP
#include <omp.h>
#endif

RVE::RVE(const std::string& config_file) {
    std::ifstream file(config_file);
    
    file >> rot_mode;
    
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
            grains.emplace_back(mat, random_q, rot_mode);
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
    
    int n = static_cast<int>(grains.size());
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n; i++) {
        grains[i].Simulate_Step(dt, L);
    }
    
    for (int i = 0; i < n; i++) {
        double weight = grains[i].get_volume() / total_volume;
        stress_av += grains[i].get_stress_lab() * weight;
    }
    
    time += dt;
    
    F_av = (dt * L).exp() * F_av;
}

void RVE::Simulate_Step_Soft(double dt, const Eigen::Matrix3d& dS_dt) {
    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using Matrix6d = Eigen::Matrix<double, 6, 6>;
    
    auto tensorToVoigt = [](const Eigen::Matrix3d& S) -> Vector6d {
        Vector6d v;
        v << S(0,0), S(1,1), S(2,2), S(1,2)+S(2,1), S(0,2)+S(2,0), S(0,1)+S(1,0);
        return v;
    };
    
    auto voigtToTensor = [](const Vector6d& v) -> Eigen::Matrix3d {
        Eigen::Matrix3d S;
        S << v(0), v(5)/2, v(4)/2,
             v(5)/2, v(1), v(3)/2,
             v(4)/2, v(3)/2, v(2);
        return S;
    };
    
    Vector6d dS_dt_v = tensorToVoigt(dS_dt);
    Vector6d stress_target_v = stress_av_old_v + dS_dt_v * dt;
    
    Vector6d Dp_v = Dp_av_old_v;
    
    int n = static_cast<int>(grains.size());
    
    for (int iter = 0; iter < max_iter; iter++) {
        // 1. Вычислить среднюю упругость (в ЛСК, Фойгт)
        Matrix6d P_av = Matrix6d::Zero();
        
#ifdef _OPENMP
        #pragma omp parallel
        {
            Matrix6d P_local = Matrix6d::Zero();
            #pragma omp for schedule(static)
            for (int i = 0; i < n; i++) {
                double w = grains[i].get_volume() / total_volume;
                P_local += w * grains[i].get_elasticity_lab();
            }
            #pragma omp critical
            P_av += P_local;
        }
#else
        for (int i = 0; i < n; i++) {
            double w = grains[i].get_volume() / total_volume;
            P_av += w * grains[i].get_elasticity_lab();
        }
#endif
        
        Matrix6d P_inv = P_av.inverse();
        
        // 2. Вычислить D
        Vector6d D_v = P_inv * dS_dt_v + Dp_v;
        
        // 3. Преобразовать в L (спин с прошлого шага)
        Eigen::Matrix3d D = voigtToTensor(D_v);
        Eigen::Matrix3d W = W_av_old;
        Eigen::Matrix3d L = D + W;
        
        // 4. Параллельный расчёт зёрен
#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (int i = 0; i < n; i++) {
            grains[i].Simulate_Step(dt, L);
        }
        
        // 5. Осреднить результаты
        Vector6d stress_v = Vector6d::Zero();
        Dp_v = Vector6d::Zero();
        Eigen::Matrix3d W_sum = Eigen::Matrix3d::Zero();
        
        for (int i = 0; i < n; i++) {
            double w = grains[i].get_volume() / total_volume;
            stress_v += w * tensorToVoigt(grains[i].get_stress_lab());
            Dp_v += w * grains[i].get_Dp_v();
            W_sum += w * 0.5 * (L - L.transpose());
        }
        
        // 6. Проверка сходимости
        double rel_error = (stress_v - stress_target_v).norm() 
                         / std::max(stress_target_v.norm(), 1e-10);
        
        if (rel_error < tol) {
            stress_av = voigtToTensor(stress_v);
            stress_av_old_v = stress_v;
            Dp_av_old_v = Dp_v;
            W_av_old = W_sum;
            W_av = W_sum;
            Dp_av = voigtToTensor(Dp_v);
            time += dt;
            F_av = (dt * L).exp() * F_av;
            return;
        }
    }
    
    std::cout << "Warning: soft loading did not converge in " << max_iter << " iterations" << std::endl;
}
