#pragma once
#include "../cp_model/Grain.hpp"
#include "cp_umat_types.hpp"

class CPUMAT {
public:
    static constexpr int NSTATV = 70;

    void init_from_props(const double* PROPS, int NPROPS);
    void load_state(const double* STATEV);
    void save_state(double* STATEV) const;
    void step(const double* DSTRAN, double DTIME, double* STRESS, double* DDSDDE, double& PNEWDT);

private:
    bool initialized_ = false;
    MaterialParams* raw_params_ = nullptr;  // heap-allocated, lifetime = singleton
    Grain* grain_ = nullptr;
    Vector6d D_v_buf_;
    Matrix6d T_buf_;
};
