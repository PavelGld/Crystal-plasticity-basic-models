#include "RVE.h"
#include <iostream>
#include <fstream>
#include <unsupported/Eigen/MatrixFunctions>

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
}

void RVE::Simulate_Step(double dt, const Eigen::Matrix3d& L) {
    stress_av.setZero();
    W_av.setZero();
    Dp_av.setZero();
    
    for (auto& grain : grains) {
        grain.Simulate_Step(dt, L);
    }
    
    double n = static_cast<double>(grains.size());
    for (const auto& grain : grains) {
        stress_av += grain.get_stress_lab() / n;
    }
    
    time += dt;
    
    F_av = (dt * L).exp() * F_av;
}
