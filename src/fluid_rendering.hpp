#pragma once
#include <vector>
#include <numeric>
#include <iostream>
#include <SFML/Graphics.hpp>
#include "simulation_state.hpp"
#include "utils.hpp"

struct FluidRenderer {
    FluidState &fluid_attributes;
    int n;
    float radius;
    float invSpacing;
    int32_t renderPattern = 0;

    float diffusionRatio = 0.75f;

    sf::RenderWindow &window;

    sf::VertexArray va{sf::PrimitiveType::Quads};
    sf::Texture circle_texture;
    sf::RenderStates fluidStates;
    std::vector<sf::Vertex> vaCopy;
    sf::Vector2f circle_texture_size;

    sf::RectangleShape cellDrawer;
    sf::CircleShape circleDrawer;

    sf::VertexArray cellVa{sf::PrimitiveType::Quads};
    sf::Texture cellTexture;
    sf::RenderStates cellStates;
    sf::Vector2f cell_texture_size;

    std::vector<float> particleColors;
    std::vector<int32_t> collisions;

    std::array<std::array<int, 3>, 100> velGradient;
    std::array<std::array<int, 3>, 4> velColorMap {{{50, 0, 255}, {225, 0, 225}, {255, 225, 100}, {255, 255, 125}}};

    std::array<std::array<int, 3>, 100> vortGradient;
    std::array<std::array<int, 3>, 4> vortColorMap {{{50, 0, 255}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}}};

    std::array<std::array<int, 3>, 100> tempgradient;
    std::array<std::array<int, 3>, 4> tempMap{{{0, 0, 0}, {204, 51, 0}, {255, 102, 0}, {255, 255, 102}}};
        // some nice gradients to put into these colorMaps:
        // velocity:
            // sunset: {50, 0, 255}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}
            // scientific: {0, 150, 255}, {0, 255, 0}, {255, 255, 0}, {255, 0, 0}
            // night ocean: {0, 0, 100},{0, 180, 255}, {100, 255, 255}, {255, 255, 255}
            // ocean: {0, 60, 130}, {0, 180, 255}, {180, 255, 255}, {240, 240, 240}
            // bright ocean: {100, 150, 255}, {100, 255, 255}, {255, 255, 255}, {255, 255, 255}
            // plasma: {30, 0, 70}, {120, 0, 255}, {255, 150, 255}, {255, 255, 255}
            // plasma neon: {0, 0, 100}, {128, 0, 255}, {255, 0, 128}, {255, 255, 255}
            // bright plasma neon: {80, 0, 225}, {128, 0, 255}, {255, 0, 128}, {255, 255, 255}
            // deep ocean: {0, 10, 60}, {0, 80, 180}, {0, 180, 220}, {200, 255, 255}
        // vorticity:
            // sunset: {50, 0, 255}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}
        // temperature:
            // fire: {0, 0, 0}, {204, 51, 0}, {255, 102, 0}, {255, 255, 102}
            // sunset: {50, 0, 255}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}

    std::vector<bool> debug_condition;

    FluidRenderer(FluidState &fas, sf::RenderWindow &w): fluid_attributes(fas), window(w) {
        int numParticles = fluid_attributes.num_particles;
        n = fluid_attributes.numY;
        radius = fluid_attributes.radius + 0.5f;
        invSpacing = 1.f / fluid_attributes.cellSpacing;

        this->collisions.resize(numParticles);

        this->va.resize(numParticles * 4);
        this->vaCopy.resize(numParticles * 4);
        circle_texture.loadFromFile("white_circle.png");
        circle_texture.generateMipmap();
        circle_texture_size =static_cast<sf::Vector2f>(circle_texture.getSize());
        for (int index = 0; index < numParticles; ++index) {
            int i = 4 * index;
            va[i].texCoords = {0.f, 0.f};
            va[i + 1].texCoords = {circle_texture_size.x, 0.f};
            va[i + 2].texCoords = {circle_texture_size.x, circle_texture_size.y};
            va[i + 3].texCoords = {0.f, circle_texture_size.y};

            vaCopy[i].texCoords = {0.f, 0.f};
            vaCopy[i + 1].texCoords = {circle_texture_size.x, 0.f};
            vaCopy[i + 2].texCoords = {circle_texture_size.x, circle_texture_size.y};
            vaCopy[i + 3].texCoords = {0.f, circle_texture_size.y};
        }
        fluidStates.texture = &circle_texture;


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

        debug_condition.resize(fluid_attributes.num_particles);
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

            int vort = static_cast<int>(fluid_attributes.calcVorticity(cellX, cellY)); // * 5 for curl

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

            int tempidx = std::min(tempgradient.size() - 1, static_cast<unsigned long long>(temp));

            color = sf::Color(tempgradient[tempidx][0], tempgradient[tempidx][1], tempgradient[tempidx][2], 255);

            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void updateVertexArrayDebug(int startIndex, int endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            int i = 4 * index;
            const float px = fluid_attributes.positions[2 * index];
            const float py = fluid_attributes.positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            sf::Color color = sf::Color::Red;

            if (debug_condition[index]) {
                color = sf::Color::Green;
            }


            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void updateVaDebugMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->updateVertexArrayDebug(start, end);
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
        window.draw(va, fluidStates);
    }

    void updateDivergenceVa(int start, int end, float maxDiv) {
        for (int index = start; index < end; ++index) {
            int i = 4 * index;
            cellVa[i].texCoords = {0.f, 0.f};
            cellVa[i + 1].texCoords = {cell_texture_size.x, 0.f};
            cellVa[i + 2].texCoords = {cell_texture_size.x, cell_texture_size.y};
            cellVa[i + 3].texCoords = {0.f, cell_texture_size.y};

            int cellID = fluid_attributes.fluid_cells[index];

            sf::Vector2f cell_pos = fluid_attributes.gridCellToPos(cellID);

            float px = cell_pos.x;
            float py = cell_pos.y;

            cellVa[i].position = {px - fluid_attributes.halfSpacing, py - fluid_attributes.halfSpacing};
            cellVa[i + 1].position = {px + fluid_attributes.halfSpacing, py - fluid_attributes.halfSpacing};
            cellVa[i + 2].position = {px + fluid_attributes.halfSpacing, py + fluid_attributes.halfSpacing};
            cellVa[i + 3].position = {px - fluid_attributes.halfSpacing, py + fluid_attributes.halfSpacing};

            float div = fabsf(fluid_attributes.u[cellID + n] - fluid_attributes.u[cellID] + fluid_attributes.v[cellID + 1] - fluid_attributes.v[cellID]);

            div /= maxDiv;

            float redScale = 255 * (div);
            float greenScale = 255 * (1.f - div);

            sf::Color color = sf::Color(redScale, greenScale, 0);

            //color = sf::Color(255 * (div > 0.5f), 255 * (div < 0.5f), 0);

            cellVa[i].color = color;
            cellVa[i + 1].color = color;
            cellVa[i + 2].color = color;
            cellVa[i + 3].color = color;
        }
    }

    void UpdateDivergenceVaMulti() {
        cellVa.resize(4 * fluid_attributes.num_fluid_cells);

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

        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_fluid_cells, [this, maxDiv](int start, int end) {
            this->updateDivergenceVa(start, end, maxDiv);
        });
    }

    void DrawDivergences() {
        window.draw(cellVa, cellStates);
    }

    void render_fluid() {
        if (renderPattern == 0) {
            UpdateVaDiffusionMulti();
        } else if (renderPattern == 1) {
            UpdateVaVelocityMulti();
        } else if (renderPattern == 2) {
            UpdateVaVorticityMulti();
        } else if (renderPattern == 3) {
            UpdateVaTemperatureMulti();
        } /*else if (renderPattern == 5) {
            updateVaDebugMulti();
        }*/

        if (renderPattern == 4) {
            UpdateDivergenceVaMulti();
            DrawDivergences();
            //drawActiveUVNodes(window);
        } else {
            DrawParticles();
        }

        //std::fill(begin(debug_condition), end(debug_condition), false);
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

    void setNextRenderPattern() {
        this->renderPattern++;
        if (this->renderPattern > 4) { // 5 if wanna use debug mode
            this->renderPattern = 0;
        }
    }

};