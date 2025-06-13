#include "Composite.hpp"

void Composite::properateWave(const Wave& wave) {
    // Iterate through each material in the composite
    
    if (mediums.empty()) {
        // Handle the case where there are no materials
        return;
    }

    mediums[0].propergateForward(wave);

}