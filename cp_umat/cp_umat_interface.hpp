#pragma once
#include "../cp_model/Grain.hpp"
#include "cp_umat_types.hpp"

class CPUMAT {
public:
    static constexpr int NSTATV = 70;

    static void init_params(const double* PROPS, int NPROPS);
    static Grain* get_grain_for_thread(const MaterialParams* params);
    static void load_state(Grain* g, const double* STATEV);
    static void save_state(const Grain* g, double* STATEV);
    static void step(Grain* g, const double* DSTRAN, double DTIME, double* STRESS, double* DDSDDE, double& PNEWDT, double& SSE, double& SPD);
};
