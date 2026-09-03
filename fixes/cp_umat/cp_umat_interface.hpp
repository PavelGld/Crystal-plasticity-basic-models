#pragma once
#include "../cp_model/Grain.hpp"
#include "cp_umat_types.hpp"

class CPUMAT {
public:
    static constexpr int NSTATV = 70;

    void init_params(const double* PROPS, int NPROPS);
    Grain* grain();
    static void load_state(Grain* g, const double* STATEV);
    static void save_state(const Grain* g, double* STATEV);
    void step(Grain* g, const double* DSTRAN, double DTIME, double* STRESS, double* DDSDDE, double& PNEWDT);

private:
    const MaterialParams* params_ = nullptr;
    Grain* grain_ = nullptr;
};