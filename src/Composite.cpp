#include "Composite.hpp"

#include <iostream>

void Composite::properateWave( Wave& wave) {
    // Iterate through each material in the composite
    
    std::cout << "Composite::properateWave() called" << std::endl;

    if (mediums.empty()) {
        // Handle the case where there are no materials
        return;
    }

    // Note zero is the top of the composite.
    mediums[0].initialWave(wave);

}

void Composite::addMaterial( Medium& material) {
    Composite::mediums.push_back(material);
}

void Composite::makeComposite( std::vector<float>& thickness,  std::vector<float>& density,
                    std::vector<float>& cp,  std::vector<float>& cs,
                    std::vector<float>& att_p,  std::vector<float>& att_s)
{
    std::cout << "Composite::makeComposite() called" << std::endl;

    if (thickness.size() != density.size() || thickness.size() != cp.size() ||
        thickness.size() != cs.size() || thickness.size() != att_p.size() ||
        thickness.size() != att_s.size()) {
        throw std::invalid_argument("All input vectors must have the same size.");
    }

    for (size_t i = 0; i < thickness.size(); ++i) {
        Medium medium(thickness[i], density[i], cp[i], cs[i], att_p[i], att_s[i]);
        addMaterial(medium);
    }
    std::cout << "Composite::makeComposite() completed with " << mediums.size() << " materials." << std::endl;
}


std::vector<Wave> Composite::getRp()
{
    std::cout << "Composite::getRp() called" << std::endl;
    // Return an empty vector or your actual data
    return {};
}

std::vector<Wave> Composite::getRs()
{
    std::cout << "Composite::getRs() called" << std::endl;
    return {};
}

std::vector<Wave> Composite::getTp()
{
    std::cout << "Composite::getTp() called" << std::endl;
    return {};
}

std::vector<Wave> Composite::getTs()
{
    std::cout << "Composite::getTs() called" << std::endl;
    return {};
}