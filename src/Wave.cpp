#include "Wave.hpp"

float Wave::omega = 0.0f;

bool Wave::operator==( Wave& other)  {
    if (type != other.type) {
        return false; // Different types cannot be equal
    }
    if (abs(angle - other.angle) >= 1e-6f) {
        return false; // Angles are not equal within a small tolerance
    }
    return true;
}

void Wave::clear() {
    p = std::complex<float>(0.0f,0.0f);
}

void Wave::accumulate( Wave& other) {
    p += other.p;
}