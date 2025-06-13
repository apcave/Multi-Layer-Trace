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
        virtual void stepWaves() = 0;

    public:
        enum class Type
        {
            Top,
            Intermediate,
            Bottom
        };

    protected:
        void matchWave(Wave& wave, std::vector<Wave>& waveVector);
        virtual void addWavesTop(Wave& tp, Wave& ts);
        virtual void addWavesBottom(Wave& tp, Wave& ts);

    private:
        Type mediumType; // Type of medium (Top, Intermediate, Bottom)

    protected:
        float t;   // Thickness Meters
        float rho; // Density
        float cp_r; // Compression Wave Speed
        float cs_r; // Shear Wave Speed
        float att_p; // Attenuation for P-wave
        float att_s; // Attenuation for S-wave

        std::complex k_p; // Complex wave number for P-wave
        std::complex k_s; // Complex wave number for S-wave

        std::vector<Wave> surface_rp; // P-wave on the surface Top or Bottom
        std::vector<Wave> surface_rs; // S-wave on the surface Top or Bottom

        // Note that mediums many only have on interface set.
        Interface boundryTop;
        Interface boundryBottom;
};