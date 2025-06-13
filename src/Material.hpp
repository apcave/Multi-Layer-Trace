#pragma once
#include "Wave.hpp"
#include "Medium.hpp"
#include "Interface.hpp"

#include <vector>
#include <complex>

class Material : public Medium
{
    public:
        Material(float thickness, float density, float cp, float cs, float att_p, float att_s,
                 MediumType type)
            : Medium(thickness, density, cp, cs, att_p, att_s), mediumType(type), mediumAbove(above), mediumBelow(below) {}

        void setAdjacentMediums(Material& above, Material& below) {
            mediumAbove = above;
            mediumBelow = below;
        }

        virtual ~Material() = default;


    protected:
        void stepWaves() : override;

        // Propagate the wave forward through the material
        void propergateForward(const Wave& wave);

        // Propagate the wave backward through the material
        void propergateBackward(const Wave& wave);

        // Add waves at the top interface
        void addWavesTop(const Wave& tp, const Wave& ts): override;

        // Add waves at the bottom interface
        void addWavesBottom(const Wave& tp, const Wave& ts): override;

    protected:
        void propergateForward(const Wave& wave);
        void propergateBackward(const Wave& wave);

    private:
        Wave attenuateForTransmission(Wave& wave_in);    
    
    private:
        float minP = 1.0e-6; // Minimum amplitude for wave propagation


        // Top side iterface waves.
        Material& mediumAbove; // Pointer to the medium above this one
        std::vector<Wave> rp_i1_f; // Reflected P-wave forward
        std::vector<Wave> rs_i1_f; // Reflected S-wave forward


        // Bottom side interface waves.
        Material& mediumBelow; // Pointer to the medium below this one
        std::vector<Wave> rp_i2_b; // Reflected P-wave backward
        std::vector<Wave> rs_i2_b; // Reflected S-wave backward
}   