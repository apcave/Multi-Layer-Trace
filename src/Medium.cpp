#include "Medium.hpp"
#include <cmath>

float Medium::eta = 1/(40.0f * M_PI * std::log10(std::exp(1.0f)));
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
    if (wave.isSurfaceWave()) {
        // If the wave is a surface wave, it should not be added to the wave vector.
        // Surface waves are handled separately.
        std::cout << "Surface wave detected, not adding to wave vector." << std::endl;
        return;
    }

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
    k_p = std::complex<float>(omega / cp, att_p / 8.68589f);

    if (isSolid ) {
        //k_s = (omega / cs ) * ( j *( att_s / eta ) + 1.0f);
        k_s = std::complex<float>(omega / cp, att_s / 8.68589f);
    }
    std::cout << "k_p : " << k_p << ", k_s: " << k_s << std::endl;
    std::cout << "eta :" << eta << std::endl;
}

void Medium::clearWaves()
{
    surface_pc.clear();
    surface_ps.clear();
}

void Medium::printWaves()
{
    std::cout << "Surface P-waves:" << std::endl;
    for (const auto& wave : surface_pc) {
        wave.print();
    }

    std::cout << "Surface S-waves:" << std::endl;
    for (const auto& wave : surface_ps) {
        wave.print();
    }
}