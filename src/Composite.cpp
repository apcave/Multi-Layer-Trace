#include "Composite.hpp"
#include "Material.hpp"

#include <iostream>

void Composite::properateWave( Wave& wave) {
    // Iterate through each material in the composite
    
    std::cout << "Composite::properateWave() called" << std::endl;

    if (mediums.size() < 2) {
        // Handle the case where there are no materials
        return;
    }

    for (auto& medium : mediums) {
        // Step the waves in the medium
        medium->calculateWaveNumbers(wave.omega);
        medium->clearWaves();
    }

    // Note zero is the top of the composite.
    mediums[0]->initialWave(wave);

    Medium* topMedium = mediums[0];
    Medium* bottomMedium = mediums[mediums.size() - 1];
    
    bool allDone = true;
    while (true) {
    //for (int i = 0; i < 2; ++i) {
        std::cout << "Stepping waves in all mediums..." << std::endl;
        for (int i = 0; i < mediums.size(); ++i) {
            // Step the waves in the medium
            Medium* medium = mediums[i];
            bool mediumDone = medium->stepWaves();
            allDone = allDone && mediumDone;
        }

        std::cout << "*************************************************" << std::endl;
        std::cout << "Reflected Waves..." << std::endl;
        topMedium->printWaves();
        std::cout << "Transmitted Waves..." << std::endl;
        bottomMedium->printWaves();
        std::cout << "*************************************************" << std::endl;
 
        if(allDone) {
            std::cout << "All mediums have completed stepping waves." << std::endl;
            break;
        }
        allDone = true; 
 
    }

    std::cout << "Composite::properateWave() completed." << std::endl;

}

void Composite::addMaterial( Medium* material) {
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
        if (i == 0 || i == thickness.size() - 1) {
            // Set the type of the first and last medium
            if (i == 0) {
                std::cout << "Creating top medium." << std::endl;
            } else {
                std::cout << "Creating bottom medium." << std::endl;
            }
            auto medium = new Medium(density[i], cp[i], cs[i], att_p[i], att_s[i]);
            addMaterial(medium);            
        } else {
            std::cout << "Creating intermediate material at index " << i << "." << std::endl;
            auto material = new Material(thickness[i], density[i], cp[i], cs[i], att_p[i], att_s[i]);
            addMaterial(material);            
        }
    }

    // Note the bottom layer bottom interface is not set or required.
    for (size_t i = 0; i < thickness.size()-1; ++i) {
        if ( i == 0 ) {
            Medium* currentMedium = mediums[i];
            Medium* belowMedium = mediums[i + 1];
            currentMedium->setLayerBelow(belowMedium);
        } else {
            Medium* currentMedium = mediums[i];
            Medium* aboveMedium = mediums[i - 1];
            Medium* belowMedium = mediums[i + 1];
            currentMedium->setLayerAbove(aboveMedium);
            currentMedium->setLayerBelow(belowMedium);
        }
    }

    std::cout << "Composite::makeComposite() completed with " << mediums.size() << " materials." << std::endl;
}


std::vector<Wave> Composite::getRp()
{
    // std::cout << "Composite::getRp() called" << std::endl;

    if (mediums.size() < 2) {
        // Handle the case where there are no materials
        return {};
    }

    Medium* topMedium = mediums[0];
    return topMedium->surface_pc;
}

std::vector<Wave> Composite::getRs()
{
    // std::cout << "Composite::getRs() called" << std::endl;

    if (mediums.size() < 2) {
        // Handle the case where there are no materials
        return {};
    }

    Medium* topMedium = mediums[0];
    return topMedium->surface_ps;
}

std::vector<Wave> Composite::getTp()
{
    // std::cout << "Composite::getTp() called" << std::endl;

    if (mediums.size() < 2) {
        // Handle the case where there are no materials
        return {};
    }

    Medium* bottomMedium = mediums[mediums.size() - 1];
    return bottomMedium->surface_pc;
}

std::vector<Wave> Composite::getTs()
{
    // std::cout << "Composite::getTs() called" << std::endl;

    if (mediums.size() < 2) {
        // Handle the case where there are no materials
        return {};
    }

    Medium* bottomMedium = mediums[mediums.size() - 1];
    return bottomMedium->surface_ps;
}