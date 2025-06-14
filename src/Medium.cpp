#include "Medium.hpp"

Medium::Medium(float thickness, float density, float cp, float cs, float att_p, float att_s)
{
    this->t = thickness;   // Thickness in meters
    this->rho = density; // Density in kg/m^3
    this->cp_r = cp; // Compression wave speed in m/s
    this->cs_r = cs; // Shear wave speed in m/s
    this->att_p = att_p; // Attenuation for P-wave in dB/m
    this->att_s = att_s; // Attenuation for S-wave in dB/m
}



void Medium::matchWave(Wave& wave, std::vector<Wave>& waveVector)
{
    // Add waves together if the have the same reflection angle.
    for (auto& w : waveVector) {
        if( wave == w ) {
            w.accumulate(wave);
            return;
        }
    }
    waveVector.push_back(wave);
}

void Medium::addWavesTop(Wave& tp, Wave& ts)
{
    matchWave(tp, surface_rp);
    matchWave(ts, surface_rs);
}


void Medium::addWavesBottom(Wave& tp, Wave& ts)
{
    matchWave(tp, surface_rp);
    matchWave(ts, surface_rs);
}

void  Medium::initialWave(Wave& wave)
{
    std::vector<Wave> ir = boundaryBottom.getSplitWaves(wave);

    addWavesBottom(ir[0], ir[1]);


}