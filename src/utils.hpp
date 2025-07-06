#pragma once

float clamp(float x, float min, float max) {
    return (x < min) * min + (x > max) * max + (x >= min && x <= max) * x;
}

template <typename T>
T sign(T x) {
    return (x < 0) * -1.f + (x > 0) * 1.f;
}

template <typename T>
int find(std::vector<T> *arr, T find) {
    size_t len = arr->size();
    for (int i = 0; i < len; ++i) {
        if ((*arr)[i] == find) {
            return i;
        }
    }
    return -1;
} 

void addValueToAverage(float& average, float newValue, int steps) {
    average = (average * (steps - 1) + newValue) / steps;
}

float Eps = 1e-6;

bool ltEpsPlus(float a, float b) {
    return a < b + Eps;
}

bool ltEpsMinus(float a, float b) {
    return a < b - Eps;
}

bool gtEpsPlus(float a, float b) {
    return a > b + Eps;
}

bool gtEpsMinus(float a, float b) {
    return a > b - Eps;
}

bool lteEpsPlus(float a, float b) {
    return a <= b + Eps;
}

bool lteEpsMinus(float a, float b) {
    return a <= b - Eps;
}

bool gteEpsPlus(float a, float b) {
    return a >= b + Eps;
}

bool gteEpsMinus(float a, float b) {
    return a >= b - Eps;
}

struct Vector2vu {
    std::array<float, 2> vu;
};