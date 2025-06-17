#include "Material.hpp"


void Material::propergateForward(Wave& wave)
{
    // The wave is inside the material at the top interface.

    // Attenuate the wave for transmission through the material.
    Wave wave_t = attenuateForTransmission(wave);

    // Split the wave at the bottom interface.
    std::vector<Wave> ir = boundaryTop.getSplitWaves(wave_t);

    addWavesBottom(ir[0], ir[1]);

    // The bottom of this material is the top of the medium below.
    mediumBelow->addWavesTop(ir[2], ir[3]); 
}



void Material::propergateBackward(Wave& wave)
{
    // The wave is inside the material at the bottom interface.

    // Attenuate the wave for transmission through the material.
    Wave wave_t = attenuateForTransmission(wave);

    // Split the wave at the top interface.
    std::vector<Wave> ir = boundaryTop.getSplitWaves(wave_t);
    
    addWavesTop(ir[0], ir[1]);

    // The top of this material is the bottom of the medium above.
    mediumAbove->addWavesBottom(ir[2], ir[3]);
}



void Material::addWavesTop(Wave& tp,  Wave& ts)
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

void Material::printWaves()
{
    float maxPressure = 0.0f;

    std::cout << "Top side interface waves:" << std::endl;
    for (const auto& wave : rp_i1_f) {
        if (abs(wave.p) > maxPressure) {
            maxPressure = abs(wave.p);
        }
        wave.print();
    }
    for (const auto& wave : rs_i1_f) {
        if (abs(wave.p) > maxPressure) {
            maxPressure = abs(wave.p);
        }
        wave.print();
    }

    std::cout << "Bottom side interface waves:" << std::endl;
    for (const auto& wave : rp_i2_b) {
        if (abs(wave.p) > maxPressure) {
            maxPressure = abs(wave.p);
        }
        wave.print();
    }
    for (const auto& wave : rs_i2_b) {
        if (abs(wave.p) > maxPressure) {
            maxPressure = abs(wave.p);
        }
        wave.print();
    }
    std::cout << "Max Pressure: " << maxPressure << std::endl;    
}

bool Material::stepWaves() 
{
    // Steps each wave between the top and bottom interfaces of the material.
    // Waves attenuate as they propagate through the material.
   // std::cout << "Material::stepWaves() called" << std::endl;

    bool allDone = true;
    for (auto& wave : rp_i1_f) {
        if (abs(wave.p) < minP) {
            std::cout << "Skipping P-wave with low pressure: " << abs(wave.p) << std::endl;
            continue;
        }

        allDone = false;
        Wave w = wave;
        wave.clear();
        propergateForward(w);
    }
    for (auto& wave : rs_i1_f) {
        if (abs(wave.p) < minP) {
            std::cout << "Skipping S-wave with low pressure: " << abs(wave.p) << std::endl;
            continue;
        }

        allDone = false;
        Wave w = wave;
        wave.clear();
        propergateForward(w);
    }

    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    printWaves();

    for (auto& wave : rp_i2_b) {
        std::cout << "p: " << abs(wave.p) << std::endl;
        if (abs(wave.p) < minP) {
            std::cout << "Skipping P-wave with low pressure: " << abs(wave.p) << std::endl;
            continue;
        }
        allDone = false;
        Wave w = wave;
        wave.clear();
        propergateBackward(w);
    }

    for (auto& wave : rs_i2_b) {
        if (abs(wave.p) < minP) {
            std::cout << "Skipping S-wave with low pressure: " << abs(wave.p) << std::endl;
            continue;
        }
        allDone = false;        
        Wave w = wave;
        wave.clear();
        propergateBackward(w);
    }

    if (allDone) {
        std::cout << "All waves in Material have completed stepping." << std::endl;
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        printWaves();        
    }



    return allDone;
}

Wave Material::attenuateForTransmission(Wave& wave_in)
{
    // Attenuate the wave for transmission through the material.
    // This is a simple model where the wave is attenuated by a factor of 1/e per wavelength.
    // The attenuation factor is based on the wave's type (P or S) and the material's properties.
    std::complex<float> j = std::complex<float>(0.0f, 1.0f);


    std::complex<float> distance  = thickness / cos(wave_in.angle);

    if ( abs(wave_in.p) > 10.0f )
    {
        std::cout << "Wave is amplifying error." << std::endl;
        throw std::runtime_error("Wave is amplifying error.");
    }

    //float x_offset = thickness * tan(real(wave_in.angle));

    // std::cout << "Attenuating wave for transmission through material." << std::endl;
    // std::cout << "In Wave  :" << std::endl;
    // wave_in.print();
    std::cout << "Transmission :" << std::endl;

    cn scale;
    Wave wave_out = wave_in;
    if (wave_in.type == Wave::Type::S) {
        scale = exp(j*k_s*distance);
        wave_out.p *= scale;
        std::cout << "k_s: " << k_s << std::endl;
    } else {
        scale = exp(j*k_p*distance);
        wave_out.p *= scale;
        std::cout << "k_p: " << k_p << std::endl;
        std::cout << "test: " << j*k_p*distance << std::endl; 
    }

    if(abs(scale) - 1 > 1e-6 || isnan(wave_out.p.real()) || isnan(wave_out.p.imag())) {
        std::cout << "Distance: " << distance << std::endl;
        std::cout << "p_in : " << wave_in.p << std::endl;
        std::cout << "p_out: " << wave_out.p << std::endl;
        
        if (wave_in.type == Wave::Type::S) {
            std::cout << "k_s: " << k_s << std::endl;
            std::cout << "Scale: " << abs(exp(-j*k_s*distance)) << std::endl;
        } else {
            std::cout << "k_p" << k_p << std::endl;
            std::cout << "Scale: " << abs(exp(-j*k_p*distance)) << std::endl;
        }        
        throw std::runtime_error("Attenuated wave pressure is larger than input wave pressure.");
    }

    

    std::cout << "Scale: " << scale << ", " << abs(scale) << std::endl;
    std::cout << "Distance: " << distance << std::endl;
    wave_in.print();
    wave_out.print();


    return wave_out;
}

void Material::clearWaves()
{
    rp_i1_f.clear();
    rs_i1_f.clear();
    rp_i2_b.clear();
    rs_i2_b.clear();
    Medium::clearWaves();
}