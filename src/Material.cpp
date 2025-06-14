#include "Material.hpp"


void Material::propergateForward(Wave& wave)
{
    // The wave is inside the material at the top interface.

    // Attenuate the wave for transmission through the material.
    Wave wave_t = attenuateForTransmission(wave);

    // Split the wave at the bottom interface.
    std::vector<Wave> ir = boundaryTop.getSplitWaves(wave);

    addWavesBottom(ir[0], ir[1]);

    // The bottom of this material is the top of the medium below.
    mediumBelow->addWavesTop(ir[2], ir[3]); 
}



void Material::propergateBackward(Wave& wave)
{
    // The wave is inside the material at the bottom interface.

    // Attenuate the wave for transmission through the material.
    wave = attenuateForTransmission(wave);

    // Split the wave at the top interface.
    std::vector<Wave> ir = boundaryTop.getSplitWaves(wave);
    
    addWavesTop(ir[0], ir[1]);

    // The top of this material is the bottom of the medium above.
    mediumAbove->addWavesBottom(ir[2], ir[3]);
}



void Material::addWavesTop(Wave& tp,  Wave& ts)
{
    // Waves transmitted through the interface are treated the same as reflected waves.    
    matchWave(tp, rp_i1_f);
    matchWave(ts, rs_i1_f);
}


void Material::addWavesBottom(Wave& tp, Wave& ts)
{
    // Waves transmitted through the interface are treated the same as reflected waves.
    matchWave(tp, rp_i2_b);
    matchWave(ts, rs_i2_b);
}



void Material::stepWaves() 
{
    // Steps each wave between the top and bottom interfaces of the material.
    // Waves attenuate as they propagate through the material.

    for (auto& wave : rp_i1_f) {
        if (abs(wave.p) < minP) {
            continue;
        }        
        Wave w = wave;
        wave.clear();
        propergateForward(w);
    }
    for (auto& wave : rs_i1_f) {
        if (abs(wave.p) < minP) {
            continue;
        }
        Wave w = wave;
        wave.clear();
        propergateForward(w);
    }

    for (auto& wave : rp_i2_b) {
        if (abs(wave.p) < minP) {
            continue;
        }
        Wave w = wave;
        wave.clear();
        propergateBackward(w);
    }

    for (auto& wave : rs_i2_b) {
        if (abs(wave.p) < minP) {
            continue;
        }        
        Wave w = wave;
        wave.clear();
        propergateBackward(w);
    }
}

Wave Material::attenuateForTransmission(Wave& wave_in)
{
    // Attenuate the wave for transmission through the material.
    // This is a simple model where the wave is attenuated by a factor of 1/e per wavelength.
    // The attenuation factor is based on the wave's type (P or S) and the material's properties.

    Wave wave_out = wave_in; 

    return wave_out;
}