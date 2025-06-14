#include "Wave.hpp"

bool Wave::operator==( Wave& other)  {
    return (type == other.type) && (abs(angle - other.angle) < 1e-6f) && (abs(p - other.p) < 1e-6f);
}

void Wave::clear() {
    angle = 0.0f;
}

void Wave::accumulate( Wave& other) {
    p += other.p;
}