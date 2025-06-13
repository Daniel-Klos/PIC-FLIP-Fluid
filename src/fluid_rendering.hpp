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

    float diffusionRatio = 0.85f;

    sf::RenderWindow &window;

    sf::VertexArray va{sf::PrimitiveType::Quads};
    sf::Texture texture;
    sf::RenderStates states;
    std::vector<sf::Vertex> vaCopy;
    sf::Vector2f texture_size;

    sf::RectangleShape cellDrawer;
    sf::CircleShape circleDrawer;

    sf::VertexArray cellVa{sf::PrimitiveType::Quads};
    sf::Texture cellTexture;
    sf::RenderStates cellStates;
    sf::Vector2f cell_texture_size;

    std::vector<float> particleColors;
    std::vector<int32_t> collisions;

    std::array<std::array<int, 3>, 100> velGradient;
    std::array<std::array<int, 3>, 4> velColorMap {{{255, 100, 0}, {255, 0, 100}, {100, 0, 255}, {0, 255, 255}}};

    std::array<std::array<int, 3>, 100> vortGradient;
    std::array<std::array<int, 3>, 4> vortColorMap {{{0, 0, 64}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}}};

    std::array<std::array<int, 3>, 100> tempgradient;
    std::array<std::array<int, 3>, 4> tempMap{{{0, 0, 0}, {204, 51, 0}, {255, 102, 0}, {255, 255, 102}}};
        // some nice gradients to put into these colorMaps:
        // velocity:
            // scientific: {0, 150, 255}, {0, 255, 0}, {255, 255, 0}, {255, 0, 0}
            // night ocean: {0, 0, 100},{0, 180, 255}, {100, 255, 255}, {255, 255, 255}
            // ocean: {0, 60, 130}, {0, 180, 255}, {180, 255, 255}, {240, 240, 240}
            // plasma: {30, 0, 70}, {120, 0, 255}, {255, 150, 255}, {255, 255, 255}
            // plasma neon: {0, 0, 100}, {128, 0, 255}, {255, 0, 128}, {255, 255, 255}
            // bright plasma neon: {80, 0, 225}, {128, 0, 255}, {255, 0, 128}, {255, 255, 255}
            // deep ocean: {0, 10, 60}, {0, 80, 180}, {0, 180, 220}, {200, 255, 255}
        // vorticity:
            // sunset: {0, 0, 64}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}
        // temperature:
            // fire: {0, 0, 0}, {204, 51, 0}, {255, 102, 0}, {255, 255, 102}

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


        cellTexture.loadFromFile("gray_square.png");
        cellTexture.generateMipmap();
        cell_texture_size = static_cast<sf::Vector2f>(cellTexture.getSize());
        cellStates.texture = &cellTexture;


        this->particleColors.resize(3 * numParticles);
        for (int i = 2; i < 3 * numParticles; i += 3) {
            particleColors[i] = 255;
        }


        auto size = sf::Vector2f(fluid_attributes.cellSpacing, fluid_attributes.cellSpacing);
        this->cellDrawer.setSize(size);
        this->cellDrawer.setOutlineColor(sf::Color::White);
        this->cellDrawer.setOutlineThickness(1.f);

        circleDrawer.setRadius(fluid_attributes.cellSpacing / 6);   


        // lerp between the values in colorMap to create a gradient array 
        float num_colors = velColorMap.size() - 1; // number of colors - 1
        float num_steps = 1.f * velGradient.size() / num_colors; //num_steps = 50 * key_range
        int index = 0;
        for (int i = 0; i < num_colors; ++i) {  
            for (int x = 0; x < num_steps; ++x) {
                float t = 1.f * x / num_steps;  // Interpolation factor
                // lerp for r, g, b values between colorMap[i] and colorMap [i+1]
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
        return curl * curl; // std::abs for just curl
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

            int vel = static_cast<int>((fluid_attributes.velocities[2 * index] * fluid_attributes.velocities[2 * index] + fluid_attributes.velocities[2 * index + 1] * fluid_attributes.velocities[2 * index + 1]) / 15000); // 15000
            
            vel = clamp(vel, 0, velGradientSize); 
            
            color = sf::Color(velGradient[vel][0], velGradient[vel][1], velGradient[vel][2], 255);

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

            int vort = static_cast<int>(calcVorticity(cellX, cellY)); // * 5 for curl

            vort = clamp(vort, 0, vortGradientSize);

            sf::Color color;

            color = sf::Color(vortGradient[vort][0], vortGradient[vort][1], vortGradient[vort][2], 255);

            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void updateVertexArrayDiffusion(const uint32_t startIndex, const uint32_t endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            const float s = 0.8f;

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

            float temp = std::abs(fluid_attributes.temperatures[index]);

            if (temp < 30 && fluid_attributes.fireActive) {
                temp = 0.f;
            }

            const int tempRadius = std::min(temp, fluid_attributes.radius);

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

    void updateVertexArrayCustom(int startIndex, int endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            int i = 4 * index;
            const float px = fluid_attributes.positions[2 * index];
            const float py = fluid_attributes.positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            sf::Color color = sf::Color::Green;

            if (index < fluid_attributes.num_particles / 2) {
                color = sf::Color::Red;
            }


            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void UpdateVaCustomMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->updateVertexArrayCustom(start, end);
        });
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


    /*void updateDivergenceVa(int start, int end) {
        for (int index = start; index < end; ++index) {
            int i = 4 * index;
            cellVa[i].texCoords = {0.f, 0.f};
            cellVa[i + 1].texCoords = {obstacle_texture_size.x, 0.f};
            cellVa[i + 2].texCoords = {obstacle_texture_size.x, obstacle_texture_size.y};
            cellVa[i + 3].texCoords = {0.f, obstacle_texture_size.y};

            cellVa[i].color = gray;
            cellVa[i + 1].color = gray;
            cellVa[i + 2].color = gray;
            cellVa[i + 3].color = gray;

            
        }
    }

    void UpdateDivergenceVaMulti() {
        cellVa.resize(fluid_attributes.num_fluid_cells);

        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_fluid_cells, [this](int start, int end) {
            this->updateDivergenceVa(start, end);
        });
    }*/

    void DrawDivergences(sf::RenderWindow& window) {
        float maxDiv = 0.f;
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int idx = i * n + j;
                if (fluid_attributes.cellType[idx] == fluid_attributes.FLUID_CELL) {
                    float div = fabsf(fluid_attributes.u[(i + 1) * n + j] - fluid_attributes.u[idx] + fluid_attributes.v[idx + 1] - fluid_attributes.v[idx]);

                    if (div > maxDiv) {
                        maxDiv = div;
                    }
                }
            }
        }

        if (maxDiv == 0.f) maxDiv = 1e-5f;

        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int idx = i * n + j;
                if (fluid_attributes.cellType[idx] != fluid_attributes.FLUID_CELL) continue;
                float div = fabsf(fluid_attributes.u[(i + 1) * n + j] - fluid_attributes.u[idx] + fluid_attributes.v[idx + 1] - fluid_attributes.v[idx]);
                cellDrawer.setPosition(i * fluid_attributes.cellSpacing, j * fluid_attributes.cellSpacing);

                div /= maxDiv;

                int redScale = 255 * (div);
                int greenScale = 255 * (1.f - div);

                cellDrawer.setFillColor(sf::Color(redScale, greenScale, 0));
                window.draw(cellDrawer);
            }
        }
    }

    void drawActiveUVNodes(sf::RenderWindow& window) {
        float dx = fluid_attributes.cellSpacing;

        for (int i = 0; i < fluid_attributes.numX; ++i) {
            for (int j = 0; j < fluid_attributes.numY; ++j) {
                int idx = i * n + j;

                bool activeSolidU = false;
                bool activeSolidV = false;
                bool airU = false;
                bool airV = false;

                if (i > 0 && ((fluid_attributes.cellType[idx] == fluid_attributes.FLUID_CELL && fluid_attributes.cellType[idx - n] == fluid_attributes.SOLID_CELL) || (fluid_attributes.cellType[idx - n] == fluid_attributes.FLUID_CELL && fluid_attributes.cellType[idx] == fluid_attributes.SOLID_CELL)))
                    activeSolidU = true;
                        
                if ((fluid_attributes.cellType[idx] == fluid_attributes.FLUID_CELL && fluid_attributes.cellType[idx - 1] == fluid_attributes.SOLID_CELL) || (fluid_attributes.cellType[idx - 1] == fluid_attributes.FLUID_CELL && fluid_attributes.cellType[idx] == fluid_attributes.SOLID_CELL))
                    activeSolidV = true;

                if (fluid_attributes.cellType[idx] != fluid_attributes.FLUID_CELL && fluid_attributes.cellType[idx - n] != fluid_attributes.FLUID_CELL) 
                    airU = true;

                if (fluid_attributes.cellType[idx] != fluid_attributes.FLUID_CELL && fluid_attributes.cellType[idx - 1] != fluid_attributes.FLUID_CELL)
                    airV = true;

                if ((fluid_attributes.u[idx] != 0 || activeSolidU) && !airU) {
                    float uX = i * dx - fluid_attributes.cellSpacing / 6;
                    float uY = (j + 0.5f) * dx - fluid_attributes.cellSpacing / 6;
                    circleDrawer.setPosition(uX, uY);
                    circleDrawer.setFillColor(sf::Color(0, 80, 255)); // green for u
                    window.draw(circleDrawer);
                }

                if ((fluid_attributes.v[idx] != 0 || activeSolidV) && !airV) {
                    float vX = (i + 0.5f) * dx - fluid_attributes.cellSpacing / 6;
                    float vY = j * dx - fluid_attributes.cellSpacing / 6;
                    circleDrawer.setPosition(vX, vY);
                    circleDrawer.setFillColor(sf::Color(255, 0, 0)); // red for v
                    window.draw(circleDrawer);
                }
            }
        }
    }

};