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
