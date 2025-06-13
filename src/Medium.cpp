#include "Medium.hpp"

void Medium::matchWave(const Wave& wave, std::vector<Wave>& waveVector)
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

void Medium::addWavesTop(const Wave& tp, Wave& ts)
{
    matchWave(tp, surface_rp);
    matchWave(ts, surface_rs);
}


void Medium::addWavesBottom(const Wave& tp, Wave& ts)
{
    matchWave(tp, surface_rp);
    matchWave(ts, surface_rs);
}

