#pragma once
#include "Wave.hpp"
#include "Interface.hpp"

#include <vector>
#include <complex>
#include <iostream>

class Medium
{
    public:
        Medium(float thickness, float density, float cp, float cs, float att_p, float att_s);

        virtual ~Medium() = default;

        // Pure virtual function to be implemented by derived classes
        virtual void stepWaves() {};

        void initialWave(Wave& wave);

        void setLayerAbove(Medium* medium);
        void setLayerBelow(Medium* medium);

    public:
        void matchWave(Wave& wave, std::vector<Wave>& waveVector);
        virtual void addWavesTop(Wave& tp, Wave& ts);
        virtual void addWavesBottom(Wave& tp, Wave& ts);


    protected:
        float t;   // Thickness Meters
        float rho; // Density
        float cp_r; // Compression Wave Speed
        float cs_r; // Shear Wave Speed
        float att_p; // Attenuation for P-wave
        float att_s; // Attenuation for S-wave

        std::complex<float> k_p; // Complex wave number for P-wave
        std::complex<float> k_s; // Complex wave number for S-wave

        // Note that mediums many only have on interface set.
        Interface boundaryTop;
        Interface boundaryBottom;

        Medium* mediumAbove; // Pointer to the medium above this one
        Medium* mediumBelow; // Pointer to the medium below this one

    public:
        std::vector<Wave> surface_pc; // P-wave on the surface of the medium
        std::vector<Wave> surface_ps; // S-wave on the surface of the medium
 
};