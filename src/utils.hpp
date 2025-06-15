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