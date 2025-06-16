#include "Medium.hpp"
#include <cmath>

float Medium::eta = 40.0f * M_PI * std::log10(std::exp(1.0f));
std::complex<float> Medium::j = std::complex<float>(0.0f, 1.0f); // Imaginary unit for complex calculations

Medium::Medium( float density, float cp, float cs, float att_p, float att_s)
{
    this->rho = density; // Density in kg/m^3
    this->cp = cp; // Compression wave speed in m/s
    this->cs = cs; // Shear wave speed in m/s
    this->att_p = att_p; // Attenuation for P-wave in dB/m
    this->att_s = att_s; // Attenuation for S-wave in dB/m

    isSolid = true;
    if (cs < 0.001f ) {
        isSolid = false;
    }
}

bool Medium::stepWaves() { 
    // std::cout << "Medium::stepWaves() called" << std::endl;
    return true;
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
    matchWave(tp, surface_pc);
    matchWave(ts, surface_ps);
}


void Medium::addWavesBottom(Wave& tp, Wave& ts)
{
    matchWave(tp, surface_pc);
    matchWave(ts, surface_ps);
}

void Medium::initialWave(Wave& wave)
{

    std::vector<Wave> ir = boundaryBottom.getSplitWaves(wave);


    addWavesTop(ir[0], ir[1]);
    mediumBelow->addWavesTop(ir[2], ir[3]);
}

void Medium::setLayerAbove(Medium* medium)
{
    // Create the boundary conditions at the top of this medium.
    mediumAbove = medium;
    boundaryTop.setFirstMedium(this);
    boundaryTop.setSecondMedium(medium);
}

void Medium::setLayerBelow(Medium* medium)
{
    // Create the boundary conditions at the bottom of this medium.
    mediumBelow = medium;
    boundaryBottom.setFirstMedium(this);
    boundaryBottom.setSecondMedium(medium);
}

void Medium::calculateWaveNumbers(float omega)
{
    //k_p = (omega / cp ) * ( j *( att_p / eta ) + 1.0f);
    k_p = omega / std::complex<float>(cp, att_p);

    if (isSolid ) {
        //k_s = (omega / cs ) * ( j *( att_s / eta ) + 1.0f);
        k_s = omega / std::complex<float>(cs, att_s);

    }
}

void Medium::clearWaves()
{
    surface_pc.clear();
    surface_ps.clear();
}
