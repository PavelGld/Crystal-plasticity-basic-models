#pragma once
#include "../cp_model/Grain.hpp"
#include "cp_umat_types.hpp"

class CPUMAT {
public:
    static constexpr int NSTATV = 70;

    void init_params(const double* PROPS, int NPROPS);
    Grain* grain() { return grain_; }
    void set_initial_orientation(int el_idx);
    static void load_state(Grain* g, const double* STATEV);
    static void save_state(const Grain* g, double* STATEV);
    void step(Grain* g, const double* DSTRAN, double DTIME, double* STRESS, double* DDSDDE, double& PNEWDT);

private:
    const MaterialParams* params_ = nullptr;
    int n_elems_ = 0;
    double* quat_data_ = nullptr;  // raw heap: n_elems * 4 doubles — no Eigen in static storage
    Grain* grain_ = nullptr;        // single reusable grain
};
