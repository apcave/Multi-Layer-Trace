#include "Interface.hpp"



std::vector<Wave> Interface::getSplitWaves( Wave& wave)
{
    Wave rp;
    Wave rs;
    Wave tp;
    Wave ts;



    std::vector<Wave> result;
    result.push_back(rp);
    result.push_back(rs);
    result.push_back(tp);
    result.push_back(ts);

    return result;
}


void Interface::setFirstMedium(Medium* medium)
{
    firstMedium = medium;
}

void Interface::setSecondMedium(Medium* medium)
{
    secondMedium = medium;
}