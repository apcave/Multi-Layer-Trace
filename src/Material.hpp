#pragma once
#include "Wave.hpp"
#include "Medium.hpp"
#include "Interface.hpp"

#include <vector>
#include <complex>

class Material : public Medium
{
    public:
        Material(float thickness, float density, float cp, float cs, float att_p, float att_s)
            : Medium(density, cp, cs, att_p, att_s) {this->thickness = thickness; }

        void setAdjacentMediums(Material* above, Material* below) {
            mediumAbove = above;
            mediumBelow = below;
        }

        virtual ~Material() = default;
        bool stepWaves() override;
        void clearWaves() override;

    protected:        

        // Propagate the wave forward through the material
        void propergateForward( Wave& wave);

        // Propagate the wave backward through the material
        void propergateBackward( Wave& wave);

        // Add waves at the top interface
        virtual void addWavesTop(Wave& tp, Wave& ts) override;

        // Add waves at the bottom interface
        virtual void addWavesBottom(Wave& tp, Wave& ts) override;

    private:
        Wave attenuateForTransmission(Wave& wave_in);    
        void printWaves();
    
    private:
        float minP = 1.0e-12; // Minimum amplitude for wave propagation
        float thickness; // Thickness of the material in meters


        // Top side iterface waves.
        std::vector<Wave> rp_i1_f; // Reflected P-wave forward
        std::vector<Wave> rs_i1_f; // Reflected S-wave forward


        // Bottom side interface waves.
        std::vector<Wave> rp_i2_b; // Reflected P-wave backward
        std::vector<Wave> rs_i2_b; // Reflected S-wave backward
};   