#include "Interface.hpp"



std::vector<Wave> Interface::getSplitWaves(const Wave& wave)
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