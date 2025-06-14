#pragma once

class Medium;

#include "Wave.hpp"
#include <vector>

class Interface {

    public:
        void setFirstMedium(Medium* medium);
        void setSecondMedium(Medium* medium);

        std::vector<Wave> getSplitWaves( Wave& wave);


    private:
        Medium* firstMedium;  // Pointer to the first medium
        Medium* secondMedium; // Pointer to the second medium

};