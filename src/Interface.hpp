#pragma once
#include "Medium.hpp"

#include <vector>
#include <complex>

class Interface {

    public:
        void setFirstMedium(Medium& medium) {
            firstMedium = medium;
        }
        void setSecondMedium(Medium& medium) {
            secondMedium = medium;
        }

        std::vector<Wave> getSplitWaves(const Wave& wave);


    private:
        Medium& firstMedium;  // Pointer to the first medium
        Medium& secondMedium; // Pointer to the second medium

}