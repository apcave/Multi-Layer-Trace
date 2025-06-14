#pragma once
#include <complex>


/*
    * Wave class respresents a pressure wave in a medium. 
    * It can be either a primary (P) or secondary (S) wave, characterized by its pressure amplitude, angle of incidence, and type.
    * Where the primary wave is compression wave and the secondary wave is shear wave.
    * 
    * The design decision wave made to not use pointers for the waves.
    * There will be many waves in the system and these very simple objects and the use of memory addresses
    * make for efficient use of the processor.
*/


class Wave
{
    enum class Type
    {
        P, // Primary wave
        S, // Secondary wave
    };

    public:
        Wave() : type(Type::P), p(0.0f, 0.0f), angle(0.0f) {} // Default ructor
        Wave(Type type, std::complex<float> p, float angle)
            : type(type), p(p), angle(angle) {}

        bool operator==( Wave& other) ;
        void clear();
        void accumulate( Wave& other);


    public:
        Type type;              // Wave type (P or S)
        std::complex<float> p;  // Pressure amplitude
        static float omega;     // Angular frequency in radians per second
        float angle;            // Angle of incidence in radians
};