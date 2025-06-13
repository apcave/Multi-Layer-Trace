#include "Material.hpp"


void Material::propergateForward(const Wave& wave)
{
    // The wave is inside the material at the top interface.

    // Attenuate the wave for transmission through the material.
    Wave wave = attenuateForTransmission(wave);

    // Split the wave at the bottom interface.
    wave = attenuateForTransmission(wave);

    // Split the wave at the bottom interface.
    vector<Wave> ir = boundaryTop.getSplitWaves(wave);

    addWaveBottom(ir[0], ir[1]);

    // The bottom of this material is the top of the medium below.
    mediumBelow->addWavesTop(ir[2], ir[3]); 
}



void Material::propergateBackward(const Wave& wave)
{
    // The wave is inside the material at the bottom interface.

    // Attenuate the wave for transmission through the material.
    wave = attenuateForTransmission(wave);

    // Split the wave at the top interface.
    vector<Wave> ir = boundaryTop.getSplitWaves(wave);
    
    addWaveTop(ir[0], ir[1]);

    // The top of this material is the bottom of the medium above.
    mediumAbove->addWavesBottom(ir[2], ir[3]);
}



void Material::addWavesTop(Wave& tp, const Wave& ts)
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


void Material::matchWave(const Wave& wave, vector<Wave>& waveVector)
{
    for (auto& w : waveVector) {
        if( wave == w ) {
            w.accumulate(wave);
            return;
        }
    }
    waveVector.push_back(wave);
}