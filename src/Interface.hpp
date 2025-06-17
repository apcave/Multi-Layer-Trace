#pragma once

class Medium;

#include "Wave.hpp"
#include <vector>

typedef std::complex<float> cn;

template<typename T>
T sq(T x) {
    return x * x;
}




class Interface {

    public:
        void setFirstMedium(Medium* medium);
        void setSecondMedium(Medium* medium);

        std::vector<Wave> getSplitWaves( Wave& wave);

        static float eta; // Attenuation factor
    private:
        Medium* firstMedium;  // Pointer to the first medium
        Medium* secondMedium; // Pointer to the second medium

        void p_solidToSolid();
        void p_fluidToFluid();
        void p_fluidToSolid();
        void p_solidToFluid();

        void s_solidToSolid();
        void s_solidToFluid();

        cn M_1;
        cn M_2;
        cn mu_1;
        cn mu_2;
        cn lambda_1;
        cn lambda_2;

        cn k_p;
        cn k_s;
        cn k_tp;
        cn k_ts;

        cn tpp, rpp;

        cn theta;
        cn theta_p, theta_s, theta_rp, theta_rs, theta_tp, theta_ts;
        cn rp, tp, rs, ts;
        float rho_1, rho_2; // Densities of the first and second mediums
};