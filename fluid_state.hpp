#pragma once
#include <vector>
#include <numeric>
#include <iostream>
#include <SFML/Graphics.hpp>
#include "fluid_state.hpp"
#include "general_math.hpp"

struct FluidRenderer {
    FluidState &fluid_attributes;
    int n;
    float radius;
    float invSpacing;

    float diffusionRatio = 0.7f;

    sf::RenderWindow &window;

    sf::VertexArray va{sf::PrimitiveType::Quads};
    sf::Texture texture;
    sf::RenderStates states;
    std::vector<sf::Vertex> vaCopy;
    sf::Vector2f texture_size;

    std::vector<float> particleColors;
    std::vector<int32_t> collisions;

    std::array<std::array<int, 3>, 100> velGradient;
    std::array<std::array<int, 3>, 4> velColorMap{{{0, 51, 102}, {0, 153, 204}, {102, 255, 204}, {255, 255, 255}}};

    std::array<std::array<int, 3>, 100> vortGradient;
    std::array<std::array<int, 3>, 4> vortColorMap{{{0, 0, 64}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}}};

    std::array<std::array<int, 3>, 100> tempgradient;
    std::array<std::array<int, 3>, 4> tempMap{{{0, 0, 0}, {204, 51, 0}, {255, 102, 0}, {255, 255, 102}}};
        // some nice gradients to put into these colorMaps:
        // scientific: {0, 150, 255}, {0, 255, 0}, {255, 255, 0}, {255, 0, 0}
        // night ocean: {0, 51, 102}, {0, 153, 204}, {102, 255, 204}, {255, 255, 255}
        // sunset: {0, 0, 64}, {128, 0, 128}, {255, 128, 0}, {255, 255, 0}
        // brighter sunset: {0, 0, 64}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}
        // orange to white: {102, 51, 0}, {204, 122, 0}, {255, 153, 51}, {255, 255, 255}
        // ice: {0, 102, 204}, {173, 216, 230}, {224, 255, 255}, {255, 250, 250}
        // lava: {128, 0, 0}, {255, 69, 0}, {255, 140, 0}, {255, 215, 0}
        // deep space: {0, 0, 32}, {64, 0, 128}, {128, 0, 255}, {192, 192, 255}
        // dark blue orange: {0, 0, 128}, {0, 128, 255}, {255, 128, 0}, {255, 255, 0}
        // lightning mcqueen: {255, 0, 0}, {255, 69, 0}, {255, 165, 0}, {255, 255, 0}
        // rainbow: {255, 0, 0}, {255, 255, 0}, {0, 255, 0}, {0, 200, 255} 

    FluidRenderer(FluidState &fas, sf::RenderWindow &w): fluid_attributes(fas), window(w) {
        int numParticles = fluid_attributes.num_particles;
        n = fluid_attributes.numY;
        radius = fluid_attributes.radius;
        invSpacing = 1.f / fluid_attributes.cellSpacing;

        this->collisions.resize(numParticles);

        this->va.resize(numParticles * 4);
        this->vaCopy.resize(numParticles * 4);
        texture.loadFromFile("white_circle.png");
        texture.generateMipmap();
        texture_size =static_cast<sf::Vector2f>(texture.getSize());
        for (int index = 0; index < numParticles; ++index) {
            int i = 4 * index;
            va[i].texCoords = {0.f, 0.f};
            va[i + 1].texCoords = {texture_size.x, 0.f};
            va[i + 2].texCoords = {texture_size.x, texture_size.y};
            va[i + 3].texCoords = {0.f, texture_size.y};

            vaCopy[i].texCoords = {0.f, 0.f};
            vaCopy[i + 1].texCoords = {texture_size.x, 0.f};
            vaCopy[i + 2].texCoords = {texture_size.x, texture_size.y};
            vaCopy[i + 3].texCoords = {0.f, texture_size.y};
        }
        states.texture = &texture;

        this->particleColors.resize(3 * numParticles);
        for (int i = 2; i < 3 * numParticles; i += 3) {
            particleColors[i] = 255;
        }

        // lerp between the values in colorMap to create a gradient array 
        float num_colors = velColorMap.size() - 1; // number of colors - 1
        float num_steps = 1.f * velGradient.size() / num_colors; //num_steps = 50 * key_range
        int index = 0;
        for (int i = 0; i < num_colors; ++i) {  
            for (int x = 0; x < num_steps; ++x) {
                float t = 1.f * x / num_steps;  // Interpolation factor
                // lerp for r, g, b values between colorMap[i] andcolorMap [i+1]
                int r = (int)(velColorMap[i][0] * (1 - t) + velColorMap[i + 1][0] * t);
                int g = (int)(velColorMap[i][1] * (1 - t) + velColorMap[i + 1][1] * t);
                int b = (int)(velColorMap[i][2] * (1 - t) + velColorMap[i + 1][2] * t);
                velGradient[index] = std::array<int, 3>{r, g, b};

                r = (int)(tempMap[i][0] * (1 - t) + tempMap[i + 1][0] * t);
                g = (int)(tempMap[i][1] * (1 - t) + tempMap[i + 1][1] * t);
                b = (int)(tempMap[i][2] * (1 - t) + tempMap[i + 1][2] * t);
                tempgradient[index] = std::array<int, 3>{r, g, b};

                r = (int)(vortColorMap[i][0] * (1 - t) + vortColorMap[i + 1][0] * t);
                g = (int)(vortColorMap[i][1] * (1 - t) + vortColorMap[i + 1][1] * t);
                b = (int)(vortColorMap[i][2] * (1 - t) + vortColorMap[i + 1][2] * t);
                vortGradient[index] = std::array<int, 3>{r, g, b};

                index++;
            }
        }
    }

    // move this into a fluidmath.hpp file
    float curl(int i, int j) {
        int idx = i * n + j;
        const float denom = 1.f / (2.f * fluid_attributes.cellSpacing);
        const int32_t leftType = fluid_attributes.cellType[idx - n] == 0;
        const int32_t rightType = fluid_attributes.cellType[idx + n] == 0;
        const int32_t topType = fluid_attributes.cellType[idx - 1] == 0;
        const int32_t bottomType = fluid_attributes.cellType[idx + 1] == 0;
        if (!leftType || !rightType || !topType || !bottomType) {
            return 0.f;
        }
        return ((fluid_attributes.v[(i + 1) * n + j] * bottomType - fluid_attributes.v[(i - 1) * n + j] * topType) - (fluid_attributes.u[i * n + j + 1] * rightType - fluid_attributes.u[i * n + j - 1] * leftType)) * denom;
    }

    float calcVorticity(int i, int j) {
        float curl = this->curl(i, j);
        return std::abs(curl);// * curl;
    }

    void updateVertexArrayVelocity(uint32_t startIndex, uint32_t endIndex) {
        int32_t velGradientSize = velGradient.size() - 1;
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            int i = 4 * index;
            const float px = fluid_attributes.positions[2 * index];
            const float py = fluid_attributes.positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            sf::Color color;

            int vel = (int)(fluid_attributes.velocities[2 * index] * fluid_attributes.velocities[2 * index] + fluid_attributes.velocities[2 * index + 1] * fluid_attributes.velocities[2 * index + 1]) / 15000; 
            
            color = sf::Color(velGradient[std::min(velGradientSize, static_cast<int32_t>(vel))][0], velGradient[std::min(velGradientSize, static_cast<int32_t>(vel))][1], velGradient[std::min(velGradientSize, static_cast<int32_t>(vel))][2], 255);

            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void updateVertexArrayVorticity(uint32_t startIndex, uint32_t endIndex) {
        int32_t vortGradientSize = vortGradient.size() - 1;
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            int i = 4 * index;
            const float px = fluid_attributes.positions[2 * index];
            const float py = fluid_attributes.positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            int cellX = px / fluid_attributes.cellSpacing;
            int cellY = py / fluid_attributes.cellSpacing;

            float vort = static_cast<float>(std::min(255, static_cast<int>(calcVorticity(cellX, cellY) * 5))); 
            sf::Color color;

            color = sf::Color(vortGradient[std::min(vortGradientSize, static_cast<int32_t>(vort))][0], vortGradient[std::min(vortGradientSize, static_cast<int32_t>(vort))][1], vortGradient[std::min(vortGradientSize, static_cast<int32_t>(vort))][2], 255);

            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void updateVertexArrayDiffusion(const uint32_t startIndex, const uint32_t endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            const float s = 1.f;

            int i = 4 * index;
            const float px = fluid_attributes.positions[2 * index];
            const float py = fluid_attributes.positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            particleColors[3 * index] = clamp(this->particleColors[3 * index] - s, 0, 255);
            particleColors[3 * index + 1] = clamp(this->particleColors  [3 * index + 1] - s, 0, 255);
            particleColors[3 * index + 2] = clamp(this->particleColors  [3 * index + 2] + s, 0, 255);

            const int xi = clamp(std::floor(px * invSpacing), 1, fluid_attributes.numX - 1);
            const int yi = clamp(std::floor(py * invSpacing), 1, fluid_attributes.numY - 1);
            const int cellNr = xi * n + yi;

            const float d0 = fluid_attributes.particleRestDensity;

            if (d0 > 0.f) {
                const float relDensity = this->fluid_attributes.cellDensities[cellNr] / d0;
                if (relDensity < diffusionRatio) { 
                    particleColors[3 * index] = 204;
                    particleColors[3 * index + 1] = 204;
                    particleColors[3 * index + 2] = 255;
                }
            }

            sf::Color color = sf::Color(particleColors[3 * index],  particleColors[3 * index + 1], particleColors[3 * index + 2]);

            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void updateVertexArrayTemperature(const uint32_t startIndex, const uint32_t endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            int i = 4 * index;
            const float px = fluid_attributes.positions[2 * index];
            const float py = fluid_attributes.positions[2 * index + 1];

            float temp = fluid_attributes.temperatures[index];

            if (temp < 30 && fluid_attributes.fireActive) {
                temp = 0.f;
            }

            const int tempRadius = std::min(temp, radius);

            va[i].position = {px - tempRadius, py - tempRadius};
            va[i + 1].position = {px + tempRadius, py - tempRadius};
            va[i + 2].position = {px + tempRadius, py + tempRadius};
            va[i + 3].position = {px - tempRadius, py + tempRadius};

            sf::Color color;

            color = sf::Color(tempgradient[std::min(tempgradient.size() - 1, static_cast<unsigned long long>(temp))][0], tempgradient[std::min(tempgradient.size() - 1, static_cast<unsigned long long>(temp))][1], tempgradient[std::min(tempgradient.size() - 1, static_cast<unsigned long long>(temp))][2], 255);

            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void UpdateVaVelocityMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->updateVertexArrayVelocity(start, end);
        });
    }

    void UpdateVaVorticityMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->updateVertexArrayVorticity(start, end);
        });
    }

    void UpdateVaDiffusionMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->updateVertexArrayDiffusion(start, end);
        });
    }

    void UpdateVaTemperatureMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->updateVertexArrayTemperature(start, end);
        });
    }

    void DrawParticles() {
        window.draw(va, states);
    }

};
