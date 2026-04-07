#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "RVE.h"
#define M_PI 3.14159265358979323846
// ========== Вывод ориентаций в файл ==========

void export_orientations(const RVE& rve, const std::string& filename, double strain, const std::string& loading_type) {
    std::ofstream out(filename);

    out << "# Loading: " << loading_type << ", Strain: " << strain << "\n";
    out << "# qw qx qy qz\n";

    for (size_t i = 0; i < rve.num_grains(); i++) {
        const Grain& grain = rve.get_grain(i);
        Eigen::Quaterniond q = grain.get_orientation();

        out << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << "\n";
    }

    out.close();
}

// Сохранение начальных ориентаций (для диагностики)
void save_initial_orientations(const RVE& rve, const std::string& filename) {
    std::ofstream out(filename);
    out << "# Initial random orientations\n";
    out << "# qw qx qy qz\n";
    for (size_t i = 0; i < rve.num_grains(); i++) {
        Eigen::Quaterniond q = rve.get_grain(i).get_orientation();
        out << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << "\n";
    }
    out.close();
}

// Вычисление среднего угла разориентации между двумя набора
double avg_misorientation(const std::string& file_initial, const std::string& file_final) {
    auto load_quats = [](const std::string& fn) {
        std::vector<Eigen::Quaterniond> qs;
        std::ifstream f(fn);
        std::string line;
        while (std::getline(f, line)) {
            if (line[0] == '#') continue;
            std::istringstream iss(line);
            double w, x, y, z;
            if (iss >> w >> x >> y >> z)
                qs.emplace_back(w, x, y, z);
        }
        return qs;
    };

    auto qi = load_quats(file_initial);
    auto qf = load_quats(file_final);

    if (qi.size() != qf.size() || qi.empty()) return -1.0;

    double sum_angle = 0.0;
    for (size_t i = 0; i < qi.size(); i++) {
        double dot = std::abs(qi[i].dot(qf[i]));
        dot = std::min(1.0, std::max(-1.0, dot));
        sum_angle += 2.0 * std::acos(dot) * 180.0 / M_PI;
    }
    return sum_angle / qi.size();
}

// ========== Моделирование с контролем деформации ==========

void run_loading(const std::string& config_file, const Eigen::Matrix3d& L, double target_strain,
                 double dt, const std::string& loading_type,
                 const std::vector<double>& checkpoints) {

    RVE rve(config_file);

    size_t current_checkpoint = 0;

    std::cout << "\n--- " << loading_type << " ---\n";
    std::cout << "RVE grains: " << rve.num_grains() << "\n";
    std::cout << "Target strain (equiv): " << target_strain * 100 << "%\n";
    std::cout << "Step size: dt = " << dt << "\n";

    int step_count = 0;
    while (rve.get_accumulated_strain() < target_strain) {
        rve.Simulate_Step(dt, L);
        step_count++;

        double acc = rve.get_accumulated_strain();

        // Проверяем checkpoint
        if (current_checkpoint < checkpoints.size()) {
            double cp = checkpoints[current_checkpoint];
            if (acc >= cp) {
                std::cout << "  Reached " << cp * 100 << "% equiv strain at step " << step_count
                          << " (actual: " << acc * 100 << "%)\n";

                std::string result_file = "orientations_" + loading_type + "_" + std::to_string(static_cast<int>(cp * 100)) + "pct.txt";
                export_orientations(rve, result_file, acc, loading_type);
                current_checkpoint++;
            }
        }

        // Прогресс
        if (step_count % 500 == 0) {
            std::cout << "  Step " << step_count << ": " << acc * 100 << "%\n";
        }
    }

    std::cout << "Completed: " << loading_type
              << " final_equiv_strain = " << rve.get_accumulated_strain() * 100 << "%\n";

    Eigen::Matrix3d stress = rve.get_stress();
    std::cout << "Final stress:\n" << stress / 1e6 << " MPa\n";

    // Проверка: покажем ориентации первых 5 зёрен
    std::cout << "Sample orientations (first 5 grains):\n";
    for (size_t i = 0; i < std::min((size_t)5, rve.num_grains()); i++) {
        auto q = rve.get_grain(i).get_orientation();
        std::cout << "  Grain " << i << ": q = [" << q.w() << ", " << q.x() << ", " << q.y() << ", " << q.z() << "]\n";
    }

    Eigen::Matrix3d Dp = rve.get_plastic_strain_rate();
    std::cout << "Plastic strain rate:\n" << Dp << "\n";
}

// ========== Rotation_validation ==========

int main() {
    std::cout << "=== Rotation Validation Test ===\n";

    const std::string config_file = "rve_validation.txt";
    const double target_strain = 1.0;  // 100%
    const double dt = 1e-4;

    // Контрольные точки
    std::vector<double> checkpoints = {0.3, 0.5, 1.0};  // 30%, 50%, 100%

    // ===== 1. Простой сдвиг (Simple Shear) =====
    // L = [[0, γ̇, 0], [0, 0, 0], [0, 0, 0]]
    //{
    //    std::cout << "\n### SIMPLE SHEAR ###\n";
    //    double gamma_dot = 0.1;
    //    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    //    L(0, 1) = gamma_dot;
	//
    //    run_loading(config_file, L, target_strain, dt, "simple_shear", checkpoints);
    //}
	//
    //// ===== 2. Одноосное растяжение (Uniaxial Tension) =====
    //// L = diag(ε̇, -ε̇/2, -ε̇/2) - изохорическое
    //{
    //    std::cout << "\n### UNIAXIAL TENSION ###\n";
    //    double eps_dot = 0.1;
    //    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    //    L(0, 0) = -0.5 * eps_dot;
    //    L(1, 1) = -0.5 * eps_dot;
    //    L(2, 2) = eps_dot;
	//
    //    run_loading(config_file, L, target_strain, dt, "uniaxial_tension", checkpoints);
    //}
	//
    //// ===== 3. Одноосное сжатие (Uniaxial Compression) =====
    //// L = diag(-ε̇, ε̇/2, ε̇/2)
    //{
    //    std::cout << "\n### UNIAXIAL COMPRESSION ###\n";
    //    double eps_dot = 0.1;
    //    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    //    L(0, 0) = 0.5 * eps_dot;
    //    L(1, 1) = 0.5 * eps_dot;
    //    L(2, 2) = -eps_dot;
	//
    //    run_loading(config_file, L, target_strain, dt, "uniaxial_compression", checkpoints);
    //}

    // ===== 4. Плоское растяжение (Plane Strain) =====
    // L = -e * k2⊗k2 + e * k3⊗k3 → L = diag(0, -ε̇, ε̇)
    // Растяжение вдоль Z, сжатие вдоль Y, X заблокирован
    {
        std::cout << "\n### PLANE STRAIN ###\n";
        double eps_dot = -0.1;
        Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
        L(1, 1) = -eps_dot;
        L(2, 2) = eps_dot;

        run_loading(config_file, L, target_strain, dt, "plane_strain", checkpoints);
    }

    std::cout << "\n=== All simulations completed ===\n";
    std::cout << "Output files: orientations_<loading>_<strain>.txt\n";
    std::cout << "Load into MTEX using quaternion format (w, x, y, z)\n";

    return 0;
}
