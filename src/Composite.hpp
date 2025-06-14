#pragma once
#include "Medium.hpp"
#include "Wave.hpp"

#include <vector>

class Composite {

    public:
        void addMaterial( Medium& material);

        void makeComposite( std::vector<float>& thickness,  std::vector<float>& density,
                           std::vector<float>& cp,  std::vector<float>& cs,
                           std::vector<float>& att_p,  std::vector<float>& att_s);

        // As a first step I will  do a depth first search with wave cut-off at a minimum amplitude.
        // Second step will be breadth first search with wave cut-off at a minimum amplitude.
        // Finally I will attempt to identify the power series expansion of the wave and use that to determine the wave propagation.
        void properateWave( Wave& incidentWave);

        std::vector<Wave> getRp();
        std::vector<Wave> getRs();
        std::vector<Wave> getTp();
        std::vector<Wave> getTs();

    private:
        std::vector<Medium> mediums;
};