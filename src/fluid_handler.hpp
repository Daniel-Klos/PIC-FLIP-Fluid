#pragma once
#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <array>

#include "collision_grid.hpp"
#include "pressure_solver.hpp"
#include "transfer_grid.hpp"
#include "fluid_rendering.hpp"
#include "scene_renderer.hpp"
#include "utils.hpp"

class FluidHandler {
    bool interacting = false;

    float checkseparationDist;
    float moveDist;

    float dragObjectX;
    float dragObjectY;
    float dragObjectPrevX;
    float dragObjectPrevY;
    sf::CircleShape dragObjectDrawer;
    float checkDragObjectseparationDist;
    float dragObjectXVel = 0;
    float dragObjectYVel = 0;

public:

    sf::CircleShape forceObjectDrawer;
    float checkForceObjectseparationDist;


    sf::RectangleShape generatorDrawer;

    int numRowsPerThread;
    int numMissedRows;

    const float colorDiffusionCoeff = 0.001f;

    float diffusionRatio;

    float scalingFactor;
    int32_t scaledWIDTH;
    int32_t scaledHEIGHT;

    CollisionGrid collisionGrid;

    FluidState &fluid_attributes;

    float dragObjectSimRadius;
    float forceObjectSimRadius = 250;
    float generatorObjectSimRadius = forceObjectSimRadius;

    float dragObjectRenderRadius;
    float forceObjectRenderRadius;
    float generatorObjectRenderRadius;

    bool forceObjectActive = true;
    bool dragObjectActive = false;
    bool generatorActive = false;

    std::vector<uint32_t> collisions;

    FluidRenderer &fluid_renderer;
    PressureSolver pressure_solver;
    TransferGrid transfer_grid;

    FluidHandler(FluidState& fas, FluidRenderer &fr): fluid_attributes(fas), fluid_renderer(fr), pressure_solver(fas), transfer_grid(fas) {

            this->collisions.resize(fluid_attributes.num_particles);

            this->scalingFactor = 2 * fluid_attributes.radius;

            this->scaledWIDTH = std::ceil(static_cast<float>(fluid_attributes.frame_context.WIDTH) / scalingFactor);
            this->scaledHEIGHT = std::ceil(static_cast<float>(fluid_attributes.frame_context.HEIGHT) / scalingFactor);

            collisionGrid = CollisionGrid(scaledWIDTH, scaledHEIGHT);

            this->moveDist = 2 * fluid_attributes.radius;
            this->checkseparationDist = moveDist * moveDist;

            this->dragObjectSimRadius = 150.f;
            this->dragObjectRenderRadius = dragObjectSimRadius / fluid_attributes.frame_context.zoom_amount;
            this->dragObjectX = std::floor(fluid_attributes.frame_context.WIDTH / 2);
            this->dragObjectY = 2 * dragObjectSimRadius + 10;
            this->dragObjectPrevX = dragObjectX;
            this->dragObjectPrevY = dragObjectY;
            this->dragObjectDrawer.setOrigin(dragObjectSimRadius, dragObjectSimRadius);
            this->dragObjectDrawer.setRadius(dragObjectSimRadius);
            this->dragObjectDrawer.setFillColor(sf::Color::Transparent);
            this->dragObjectDrawer.setOutlineThickness(1.f);
            this->dragObjectDrawer.setOutlineColor(sf::Color(255, 0, 0));
            this->checkDragObjectseparationDist = (fluid_attributes.radius + dragObjectSimRadius) * (fluid_attributes.radius + dragObjectSimRadius);

            this->forceObjectRenderRadius = forceObjectSimRadius / fluid_attributes.frame_context.zoom_amount;
            this->forceObjectDrawer.setOrigin(forceObjectSimRadius, forceObjectSimRadius);
            this->forceObjectDrawer.setRadius(forceObjectSimRadius);
            this->forceObjectDrawer.setOutlineThickness(1.f);
            this->forceObjectDrawer.setFillColor(sf::Color::Transparent);
            this->forceObjectDrawer.setOutlineColor(sf::Color::Red); 
            this->checkForceObjectseparationDist = (fluid_attributes.radius + forceObjectSimRadius) * (fluid_attributes.radius + forceObjectSimRadius);
            
            this->generatorObjectRenderRadius = generatorObjectSimRadius / fluid_attributes.frame_context.zoom_amount;
            this->generatorDrawer.setOrigin(generatorObjectSimRadius, generatorObjectSimRadius);
            this->generatorDrawer.setSize(sf::Vector2f(2 * generatorObjectSimRadius, 2 * generatorObjectSimRadius));
            this->generatorDrawer.setOutlineThickness(1.f);
            this->generatorDrawer.setFillColor(sf::Color::Transparent);
            this->generatorDrawer.setOutlineColor(sf::Color::Red); 
    }

    void createRandomPositions() {
        std::uniform_int_distribution<int> randWidth(fluid_attributes.radius, fluid_attributes.frame_context.WIDTH - fluid_attributes.radius);
        std::uniform_int_distribution<int> randHeight(fluid_attributes.radius, fluid_attributes.frame_context.HEIGHT - fluid_attributes.radius);

        std::random_device rd;
        std::mt19937 mt(rd());

        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            fluid_attributes.positions[i * 2] = randWidth(mt);
            fluid_attributes.positions[i * 2 + 1] = randHeight(mt);
        }
    }

    void integrate(const uint32_t startIndex, const uint32_t endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            fluid_attributes.positions[2 * i] += fluid_attributes.velocities[2 * i] * fluid_attributes.frame_context.dt;
            fluid_attributes.positions[2 * i + 1] += fluid_attributes.velocities[2 * i + 1] * fluid_attributes.frame_context.dt;
            fluid_attributes.velocities[2 * i] += fluid_attributes.gravityX * fluid_attributes.frame_context.dt;
            fluid_attributes.velocities[2 * i + 1] += fluid_attributes.gravityY * fluid_attributes.frame_context.dt;

            /*if (fluid_attributes.densities[i] > 0) {
                float buoyancy = fluid_attributes.gravityY * (1 - (fluid_attributes.particle_densities[i] / fluid_attributes.densities[i]));
                //float buoyancy = fluid_attributes.gravityY * (fluid_attributes.densities[i] - fluid_attributes.particle_densities[i]) / fluid_attributes.particleRestDensity;
                fluid_attributes.velocities[2 * i + 1] += buoyancy * dt;
            }*/
        }
    }

    void integrateMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->integrate(start, end);
        });
    }

    void makeFire(int start, int end) {
        for (int i = start; i < end; ++i) {
            if (fluid_attributes.positions[2 * i + 1] < fluid_attributes.frame_context.HEIGHT - fluid_attributes.cellSpacing - 10) {
                fluid_attributes.velocities[2 * i + 1] -= fluid_attributes.fireStrength * fluid_attributes.temperatures[i] * fluid_attributes.frame_context.dt;
                if (fluid_attributes.temperatures[i] > 0) {
                    fluid_attributes.temperatures[i] -= fluid_attributes.tempDiffusion * fluid_attributes.frame_context.dt;
                }
                if (collisions[i] == 0) {
                    fluid_attributes.temperatures[i] = 0;
                }
            }
        }
    }

    void makeFireMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->makeFire(start, end);
        });
    }

    void addObjectsToGrids()
    {

        collisionGrid.clear();

        const float minX = fluid_attributes.cellSpacing;
        const float maxX = fluid_attributes.frame_context.WIDTH - fluid_attributes.cellSpacing;
        const float minY = fluid_attributes.cellSpacing;
        const float maxY = fluid_attributes.frame_context.HEIGHT - fluid_attributes.cellSpacing;

        uint32_t i = 0;

        for (int32_t index = 0; index < fluid_attributes.num_particles; ++index) {
            float x = fluid_attributes.positions[2 * index];
            float y = fluid_attributes.positions[2 * index + 1];
            if (x > minX && x < maxX && y > minY && y < maxY) {
                int32_t gridX = x / scalingFactor;
                int32_t gridY = y / scalingFactor;
                collisionGrid.addAtom(gridX, gridY, i);
            }
            ++i;
        }
    }

    void solveContact(uint32_t index, uint32_t otherIndex)
    {
        constexpr float eps = 0.0001f;
        const float o2_o1X = fluid_attributes.positions[2 * index] - fluid_attributes.positions[2 * otherIndex];
        const float o2_o1Y = fluid_attributes.positions[2 * index + 1] - fluid_attributes.positions[2 * otherIndex + 1];

        const float dist2 = o2_o1X * o2_o1X + o2_o1Y * o2_o1Y;

        if (dist2 < checkseparationDist && dist2 > eps) {
            const float dist          = sqrt(dist2);
            const float delta = 0.5f * (moveDist - dist) / dist;
            const float col_vecX = o2_o1X * delta;
            const float col_vecY = o2_o1Y * delta;

            fluid_attributes.positions[2 * index] += col_vecX;
            fluid_attributes.positions[2 * index + 1] += col_vecY;

            fluid_attributes.positions[2 * otherIndex] -= col_vecX;
            fluid_attributes.positions[2 * otherIndex + 1] -= col_vecY;

            const float transfer = (fluid_attributes.temperatures[index] - fluid_attributes.temperatures[otherIndex]) * fluid_attributes.interConductivity * fluid_attributes.frame_context.dt;
            fluid_attributes.temperatures[index] -= transfer * fluid_attributes.frame_context.dt;
            fluid_attributes.temperatures[otherIndex] += transfer * fluid_attributes.frame_context.dt;

            collisions[index]++;
            collisions[otherIndex]++;
        }
    }

    void checkCellCollisions(uint32_t atom_idx, const CollisionCell& c)
    {
        for (uint32_t i = 0; i < c.objects_count; ++i) {
            solveContact(atom_idx, c.objects[i]);
        }
    }

    void processCell(const CollisionCell& c, uint32_t index)
    {
        for (uint32_t i = 0; i < c.objects_count; ++i) {
            const uint32_t idx = c.objects[i];
            for (int32_t side = 0; side < 2; ++side) {
                checkCellCollisions(idx, collisionGrid.data[index + collisionGrid.height + side]);
            }
            for (int32_t side = 0; side < 2; ++side) {
                checkCellCollisions(idx, collisionGrid.data[index + side]);   
            }
            //checkCellCollisions(idx, grid.data[index - grid.height]);
            checkCellCollisions(idx, collisionGrid.data[index - collisionGrid.height + 1]);
        }
    }

    void solveCollisionThreaded(uint32_t start, uint32_t end)
    {
        for (uint32_t idx = start; idx < end; ++idx) {
            processCell(collisionGrid.data[idx], idx);
        }
    }

    void makeForceObjectQueries(const int32_t strength) {
        const int32_t numCovered = std::ceil(forceObjectSimRadius / scalingFactor);

        const uint32_t mouseColumn = std::floor(fluid_attributes.frame_context.world_mouse_pos.x / scalingFactor);
        const uint32_t mouseRow = std::floor(fluid_attributes.frame_context.world_mouse_pos.y / scalingFactor);

        for (int32_t i = -numCovered; i < numCovered + 1; ++i) {
            for (int32_t j = -numCovered; j < numCovered + 1; ++j) {

                if (mouseRow + i <= 1 || mouseRow + i >= scaledHEIGHT - 1 || mouseColumn + j <= 1 || mouseColumn + j >= scaledWIDTH - 1) continue;

                const auto cell = collisionGrid.data[mouseRow + i + collisionGrid.height * (mouseColumn + j)];

                for (uint32_t id{0}; id < cell.objects_count; ++id) {
                    const uint32_t particleIndex = cell.objects[id];

                    float dx = fluid_attributes.positions[particleIndex * 2] - fluid_attributes.frame_context.world_mouse_pos.x;
                    float dy = fluid_attributes.positions[particleIndex * 2 + 1] - fluid_attributes.frame_context.world_mouse_pos.y;
                    float d2 = dx * dx + dy * dy;

                    if (d2 > checkForceObjectseparationDist || d2 == 0.0f) continue;

                    float d = std::sqrt(d2);

                    float edgeT = d / forceObjectSimRadius;
            
                    float centerT = 1 - edgeT;

                    fluid_attributes.velocities[2 * particleIndex] += (dx * strength - fluid_attributes.velocities[2 * particleIndex]) * centerT * fluid_attributes.frame_context.dt;
                    fluid_attributes.velocities[2 * particleIndex + 1] += (dy * strength - fluid_attributes.velocities[2 * particleIndex + 1]) * centerT * fluid_attributes.frame_context.dt;
                }
            }
        }
    }

    void solveCollisions()
    {
        const uint32_t slice_count  = fluid_attributes.numThreads;
        const uint32_t slice_size   = (collisionGrid.width / slice_count) * collisionGrid.height;
        const uint32_t last_cell    = fluid_attributes.numThreads * slice_size;
        
        for (uint32_t i = 0; i < fluid_attributes.numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([this, i, slice_size]{
                uint32_t const start{i * slice_size};
                uint32_t const end  {start + slice_size};
                solveCollisionThreaded(start, end);
            });
        }
        
        if (last_cell < collisionGrid.data.size()) {
            fluid_attributes.thread_pool.addTask([this, last_cell]{
                solveCollisionThreaded(last_cell, static_cast<uint32_t>(collisionGrid.data.size()));
            });
        }
        fluid_attributes.thread_pool.waitForCompletion();
    }

    void heatGround(int32_t start, int32_t end) {
        const float percentRemoved = 0.1f;
        const float leftEdge = fluid_attributes.frame_context.WIDTH * percentRemoved;
        const float rightEdge = fluid_attributes.frame_context.WIDTH * (1.0f - percentRemoved);

        for (int i = start; i < end; ++i) {
            if (fluid_attributes.positions[2 * i + 1] + fluid_attributes.radius > fluid_attributes.frame_context.HEIGHT - 2 * fluid_attributes.cellSpacing) {
                const float remove = 10.f; // % of the floor from the sides that you want not heated
                if (fluid_renderer.renderPattern == 3 && fluid_attributes.temperatures[i] < fluid_renderer.tempgradient.size() && !remove || (fluid_attributes.positions[2 * i] > leftEdge && fluid_attributes.positions[2 * i] < rightEdge)) {
                    fluid_attributes.temperatures[i] += fluid_attributes.groundConductivity * fluid_attributes.frame_context.dt;
                }
            }
        }
    }

    void heatGroundMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->heatGround(start, end);
        });
    }

    void calcVorticityConfinement(bool red, int32_t startColumn, int32_t endColumn) {
        for (int32_t i = startColumn; i < endColumn; ++i) {
            for (int32_t j = 1; j < fluid_attributes.numY - 1; ++j) {
                if (red) {
                    if ((i + j) % 2 != 0) continue; 
                }
                else {
                    if ((i + j) % 2 == 0) continue; 
                }
                float dx = abs(fluid_attributes.curl(i, j + 1)) - abs(fluid_attributes.curl(i, j - 1));
                float dy = abs(fluid_attributes.curl(i + 1, j)) - abs(fluid_attributes.curl(i - 1, j));

                const float len = std::sqrt(dx * dx + dy * dy);
                const float invLen = 1.f / (len + (len == 0.f)) - (len == 0.f);

                dx *= invLen;
                dy *= invLen;

                const float c = fluid_attributes.curl(i, j);

                fluid_attributes.v[i * fluid_attributes.n + j] += c * dx * fluid_attributes.frame_context.dt * fluid_attributes.vorticityStrength;
                fluid_attributes.u[i * fluid_attributes.n + j] += c * dy * fluid_attributes.frame_context.dt * fluid_attributes.vorticityStrength;
            }
        }
    }

    void applyVorticityConfinementRedBlack() {
        const int32_t numColumnsPerThread = (fluid_attributes.numX - 2) / fluid_attributes.numThreads;
        const int32_t numMissedColumns = fluid_attributes.numX - 2 - numColumnsPerThread * fluid_attributes.numThreads;
        
        for (int i = 0; i < fluid_attributes.numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->calcVorticityConfinement(true, 1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->calcVorticityConfinement(true, fluid_attributes.numX - 1 - numMissedColumns, fluid_attributes.numX - 1);

        fluid_attributes.thread_pool.waitForCompletion();

        for (int i = 0; i < fluid_attributes.numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->calcVorticityConfinement(false, 1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->calcVorticityConfinement(false, fluid_attributes.numX - 1 - numMissedColumns, fluid_attributes.numX - 1);

        fluid_attributes.thread_pool.waitForCompletion();
    }

    void includeDragObject() {
        // Bug when you hold force object attract then switch to drag object while holding mouse down
        const float extend = (20 * fluid_attributes.cellSpacing) / fluid_attributes.frame_context.zoom_amount;
        if (fluid_attributes.frame_context.leftMouseDown) {
            dragObjectX = fluid_attributes.frame_context.world_mouse_pos.x;
            dragObjectY = fluid_attributes.frame_context.world_mouse_pos.y;
            float vx = (dragObjectX - dragObjectPrevX) * fluid_attributes.frame_context.maxFps * !fluid_attributes.frame_context.justPressed;
            float vy = (dragObjectY - dragObjectPrevY) * fluid_attributes.frame_context.maxFps * !fluid_attributes.frame_context.justPressed;
            int objectCellX = dragObjectX / fluid_attributes.cellSpacing;
            int objectCellY = dragObjectY / fluid_attributes.cellSpacing;
            int objectCellRadius = std::ceil(dragObjectSimRadius / fluid_attributes.cellSpacing);
            for (int i = objectCellX - objectCellRadius; i < objectCellX + objectCellRadius; i++) {
                for (int j = objectCellY - objectCellRadius; j < objectCellY + objectCellRadius; j++) {
                    int cellNr = i * fluid_attributes.n + j;
                    bool notInBounds = i < 0 || i > fluid_attributes.numX || j < 0 || j > fluid_attributes.numY;
                    if (notInBounds || (fluid_attributes.cellType[i * fluid_attributes.n + j] == fluid_attributes.SOLID_CELL)) continue;
                    float dx = (i + 0.5) * fluid_attributes.cellSpacing - dragObjectX;
                    float dy = (j + 0.5) * fluid_attributes.cellSpacing - dragObjectY;

                    if (dx * dx + dy * dy < dragObjectSimRadius * dragObjectSimRadius + extend) {
                        
                        if (fluid_attributes.cellType[cellNr - fluid_attributes.n] != fluid_attributes.SOLID_CELL) {
                            fluid_attributes.u[cellNr] = vx;
                        }
                        if (fluid_attributes.cellType[cellNr + fluid_attributes.n] != fluid_attributes.SOLID_CELL) {
                            fluid_attributes.u[cellNr + fluid_attributes.n] = vx;
                        }
                        if (fluid_attributes.cellType[cellNr - 1] != fluid_attributes.SOLID_CELL) {
                            fluid_attributes.v[cellNr] = vy;
                        }
                        if (fluid_attributes.cellType[cellNr + 1] != fluid_attributes.SOLID_CELL) {
                            fluid_attributes.v[cellNr + 1] = vy;
                        }
                    }
                }
            }
            dragObjectPrevX = dragObjectX;
            dragObjectPrevY = dragObjectY;
        }
    }

    void generate() {
        float separation = fluid_attributes.radius * 2.1;
        int wideNum = std::floor((2 * generatorObjectSimRadius) / (separation));
        int highNum = wideNum;
        int numPotentiallyAdded = wideNum * highNum;

        float generatorRightBound = fluid_attributes.frame_context.world_mouse_pos.x + generatorObjectSimRadius;
        float generatorBottomBound = fluid_attributes.frame_context.world_mouse_pos.y + generatorObjectSimRadius;

        float starting_px = std::max(fluid_attributes.frame_context.world_mouse_pos.x - generatorObjectSimRadius + fluid_attributes.radius, fluid_attributes.cellSpacing + fluid_attributes.radius);
        float starting_py = std::max(fluid_attributes.frame_context.world_mouse_pos.y - generatorObjectSimRadius + fluid_attributes.radius, fluid_attributes.cellSpacing + fluid_attributes.radius);
        float px = starting_px;
        float py = starting_py;
        bool offset = true;

        fluid_attributes.positions.resize(2 * (fluid_attributes.num_particles + numPotentiallyAdded));
        int addedTo = 0;

        for (int i = 0; i < numPotentiallyAdded; ++i) {
            float prevPx = px;
            float prevPy = py;
            if (prevPy > generatorBottomBound || prevPy + fluid_attributes.radius > fluid_attributes.frame_context.HEIGHT - fluid_attributes.cellSpacing) {
                break;
            }

            int cellX = px / fluid_attributes.cellSpacing;
            int cellY = py / fluid_attributes.cellSpacing;
            int cellNr = cellX * fluid_attributes.n + cellY;

            px += separation;
            if (px > generatorRightBound) {
                px = starting_px;
                /*if (offset) {
                    px += 0.5 * separation;
                }*/
                py += separation;
                offset = !offset;
            }

            if (cellX > fluid_attributes.numX || (fluid_attributes.cellType[cellNr] != fluid_attributes.AIR_CELL || prevPx - fluid_attributes.radius < fluid_attributes.cellSpacing || prevPx + fluid_attributes.radius > fluid_attributes.frame_context.WIDTH - fluid_attributes.cellSpacing || prevPy - fluid_attributes.radius < fluid_attributes.cellSpacing)) {
                continue;
            }


            fluid_attributes.positions[2 * (fluid_attributes.num_particles + addedTo)] = prevPx;
            fluid_attributes.positions[2 * (fluid_attributes.num_particles + addedTo) + 1] = prevPy;

            addedTo++;
        }

        fluid_attributes.num_particles += addedTo;

        /*fluid_attributes.particle_densities.resize(fluid_attributes.num_particles);
        fluid_attributes.densities.resize(fluid_attributes.num_particles);*/
        fluid_attributes.positions.resize(2 * fluid_attributes.num_particles);
        fluid_attributes.velocities.resize(2 * fluid_attributes.num_particles);
        fluid_renderer.particleColors.resize(3 * fluid_attributes.num_particles);
        fluid_renderer.va.resize(4 * fluid_attributes.num_particles);

        int start = fluid_attributes.num_particles - addedTo;

        for (int i = start; i < fluid_attributes.num_particles; i++) {
            int idx1 = 2 * i;
            int idx2 = 3 * i;
            int idx3 = 4 * i;

            //fluid_attributes.particle_densities[i] = (i < fluid_attributes.num_particles / 2) ? 0.9f : 0.2f;

            fluid_attributes.velocities[idx1] = 0.f;
            fluid_attributes.velocities[idx1 + 1] = 0.f;

            fluid_renderer.particleColors[idx2] = 255;
            fluid_renderer.particleColors[idx2 + 1] = 255;
            fluid_renderer.particleColors[idx2 + 2] = 255;

            fluid_renderer.va[idx3].texCoords = {0.f, 0.f};
            fluid_renderer.va[idx3 + 1].texCoords = {fluid_renderer.circle_texture_size.x, 0.f};
            fluid_renderer.va[idx3 + 2].texCoords = {fluid_renderer.circle_texture_size.x, fluid_renderer.circle_texture_size.y};
            fluid_renderer.va[idx3 + 3].texCoords = {0.f, fluid_renderer.circle_texture_size.y};
        }

        fluid_attributes.particlesPerThread = fluid_attributes.num_particles / fluid_attributes.numThreads;
        fluid_attributes.numMissedParticles = fluid_attributes.num_particles - fluid_attributes.numThreads * fluid_attributes.particlesPerThread;

        this->collisions.resize(fluid_attributes.num_particles);
        fluid_attributes.temperatures.resize(fluid_attributes.num_particles);

        transfer_grid.topLeftCells.resize(2 * fluid_attributes.num_particles);
        transfer_grid.topRightCells.resize(2 * fluid_attributes.num_particles);
        transfer_grid.bottomLeftCells.resize(2 * fluid_attributes.num_particles);
        transfer_grid.bottomRightCells.resize(2 * fluid_attributes.num_particles);

        transfer_grid.topLeftWeights.resize(2 * fluid_attributes.num_particles);
        transfer_grid.topRightWeights.resize(2 * fluid_attributes.num_particles);
        transfer_grid.bottomLeftWeights.resize(2 * fluid_attributes.num_particles);
        transfer_grid.bottomRightWeights.resize(2 * fluid_attributes.num_particles);

        fluid_attributes.affineMats.resize(4 * fluid_attributes.num_particles); 
        fluid_attributes.dxLefts.resize(fluid_attributes.num_particles);
        fluid_attributes.dxRights.resize(fluid_attributes.num_particles);
        fluid_attributes.dyBottoms.resize(fluid_attributes.num_particles);
        fluid_attributes.dyTops.resize(fluid_attributes.num_particles);
    }

    void remove() {
        const int32_t numCovered = std::ceil(generatorObjectSimRadius / scalingFactor);
        const uint32_t mouseColumn = std::floor(fluid_attributes.frame_context.world_mouse_pos.x / scalingFactor);
        const uint32_t mouseRow = std::floor(fluid_attributes.frame_context.world_mouse_pos.y / scalingFactor);
        const size_t len = fluid_attributes.temperatures.size();
        const size_t double_len = 2 * len;
        const size_t triple_len = fluid_renderer.particleColors.size();
        const size_t quadruple_len = 2 * double_len;
    
        size_t vaSize = fluid_renderer.va.getVertexCount();
        fluid_renderer.vaCopy.resize(vaSize);
    
        for (int i = 0; i < vaSize; ++i) {
            fluid_renderer.vaCopy[i] = fluid_renderer.va[i];
        }

        size_t remove = 0;
        
        for (int32_t i = -numCovered; i < numCovered + 1; ++i) {
            for (int32_t j = -numCovered; j < numCovered + 1; ++j) {
                if (mouseRow + j <= 1 || mouseRow + j >= scaledHEIGHT - 1 || mouseColumn + i <= 1 || mouseColumn + i >= scaledWIDTH - 1)
                    continue;
    
                const auto cell = collisionGrid.data[mouseRow + j + collisionGrid.height * (mouseColumn + i)];
    
                for (uint32_t id{0}; id < cell.objects_count; ++id) {
                    const uint32_t particleIndex = cell.objects[id];

                    size_t doubleid = 2 * particleIndex;
                    size_t tripleid = 3 * particleIndex;
                    size_t quadrupleid = 4 * particleIndex;

                    size_t doubleremove = 2 * remove;
                    size_t tripleremove = 3 * remove;
                    size_t quadrupleremove = 4 * remove;
    
                    if (doubleid + 2 < double_len) {
                        /*fluid_attributes.particle_densities[particleIndex] = fluid_attributes.particle_densities[particleIndex - 1 - remove];
                        fluid_attributes.densities[particleIndex] = fluid_attributes.particle_densities[particleIndex - 1 - remove];*/

                        fluid_attributes.temperatures[particleIndex] = fluid_attributes.temperatures[len - 1 - remove];

                        fluid_attributes.positions[doubleid] = fluid_attributes.positions[double_len - 2 - doubleremove];
                        fluid_attributes.positions[doubleid + 1] = fluid_attributes.positions[double_len - 1 - doubleremove];

                        fluid_attributes.velocities[doubleid] = fluid_attributes.velocities[double_len - 2 - doubleremove];
                        fluid_attributes.velocities[doubleid + 1] = fluid_attributes.velocities[double_len - 1 - doubleremove];

                        fluid_renderer.particleColors[tripleid] = fluid_renderer.particleColors[triple_len - 3 - tripleremove];
                        fluid_renderer.particleColors[tripleid + 1] = fluid_renderer.particleColors[triple_len - 2 - tripleremove];
                        fluid_renderer.particleColors[tripleid + 2] = fluid_renderer.particleColors[triple_len - 1 - tripleremove];

                        fluid_renderer.vaCopy[quadrupleid] = fluid_renderer.vaCopy[quadruple_len - 4 - quadrupleremove];
                        fluid_renderer.vaCopy[quadrupleid + 1] = fluid_renderer.vaCopy[quadruple_len - 3 - quadrupleremove];
                        fluid_renderer.vaCopy[quadrupleid + 2] = fluid_renderer.vaCopy[quadruple_len - 2 - quadrupleremove];
                        fluid_renderer.vaCopy[quadrupleid + 3] = fluid_renderer.vaCopy[quadruple_len - 1 - quadrupleremove];
                    }

                    remove++;
                }
            }
        }

        fluid_attributes.num_particles -= remove;
    
        /*fluid_attributes.particle_densities.resize(fluid_attributes.num_particles);
        fluid_attributes.densities.resize(fluid_attributes.num_particles);*/
        fluid_attributes.positions.resize(2 * fluid_attributes.num_particles);
        fluid_attributes.velocities.resize(2 * fluid_attributes.num_particles);
        fluid_renderer.particleColors.resize(3 * fluid_attributes.num_particles);
        fluid_renderer.vaCopy.resize(4 * fluid_attributes.num_particles);
    
        fluid_renderer.va.resize(4 * fluid_attributes.num_particles);
        for (int i = 0; i < 4 * fluid_attributes.num_particles; ++i) {
            fluid_renderer.va[i] = fluid_renderer.vaCopy[i];
        }
    
        fluid_attributes.particlesPerThread = fluid_attributes.num_particles / fluid_attributes.numThreads;
        fluid_attributes.numMissedParticles = fluid_attributes.num_particles - fluid_attributes.numThreads * fluid_attributes.particlesPerThread;
    
        this->collisions.resize(fluid_attributes.num_particles);
        fluid_attributes.temperatures.resize(fluid_attributes.num_particles);
    
        transfer_grid.topLeftCells.resize(2 * fluid_attributes.num_particles);
        transfer_grid.topRightCells.resize(2 * fluid_attributes.num_particles);
        transfer_grid.bottomLeftCells.resize(2 * fluid_attributes.num_particles);
        transfer_grid.bottomRightCells.resize(2 * fluid_attributes.num_particles);
    
        transfer_grid.topLeftWeights.resize(2 * fluid_attributes.num_particles);
        transfer_grid.topRightWeights.resize(2 * fluid_attributes.num_particles);
        transfer_grid.bottomLeftWeights.resize(2 * fluid_attributes.num_particles);
        transfer_grid.bottomRightWeights.resize(2 * fluid_attributes.num_particles);
    }

    void drawDragObject(sf::RenderWindow& window) {
        this->dragObjectDrawer.setPosition(fluid_attributes.frame_context.screen_mouse_pos.x, fluid_attributes.frame_context.screen_mouse_pos.y);
        window.draw(this->dragObjectDrawer);
    }

    void drawGenerator(sf::RenderWindow& window) {
        generatorDrawer.setPosition(fluid_attributes.frame_context.screen_mouse_pos.x, fluid_attributes.frame_context.screen_mouse_pos.y);
        window.draw(generatorDrawer);
    }

    void drawForceObject(sf::RenderWindow& window) {
        forceObjectDrawer.setPosition(fluid_attributes.frame_context.screen_mouse_pos.x, fluid_attributes.frame_context.screen_mouse_pos.y);
        window.draw(forceObjectDrawer); 
    }

    void render_objects(sf::RenderWindow& window) {
        if (forceObjectActive) {
            this->drawForceObject(window);
        }
        else if (dragObjectActive) {
            this->drawDragObject(window);
        }
        else if (generatorActive) {
            this->drawGenerator(window);
        }
    }

    void addToDragObjectRadius(float add) {
        if (dragObjectRenderRadius + add > 0) {
            dragObjectRenderRadius += add;
            dragObjectDrawer.setOrigin(dragObjectRenderRadius, dragObjectRenderRadius);
            dragObjectDrawer.setRadius(dragObjectRenderRadius);

            dragObjectSimRadius = dragObjectRenderRadius / fluid_attributes.frame_context.zoom_amount;
        }
    }

    void addToForceObjectRadius(float add) {
        if (forceObjectRenderRadius + add > 0) {
            forceObjectRenderRadius += add;
            forceObjectDrawer.setOrigin(forceObjectRenderRadius, forceObjectRenderRadius);
            forceObjectDrawer.setRadius(forceObjectRenderRadius);

            float zoom = fluid_attributes.frame_context.zoom_amount;
            float zoomed_particle_rad = fluid_attributes.radius / zoom;
            forceObjectSimRadius = forceObjectRenderRadius / zoom;
            checkForceObjectseparationDist = (zoomed_particle_rad + forceObjectSimRadius) * (zoomed_particle_rad + forceObjectSimRadius);
        }
    }

    void addToGeneratorRadius(float add) {
        if (generatorObjectRenderRadius + add > 0) {
            generatorObjectRenderRadius += add;
            generatorDrawer.setOrigin(generatorObjectRenderRadius, generatorObjectRenderRadius);
            generatorDrawer.setSize(sf::Vector2f(2 * generatorObjectRenderRadius, 2 * generatorObjectRenderRadius));

            generatorObjectSimRadius = generatorObjectRenderRadius / fluid_attributes.frame_context.zoom_amount;
        }
    }

    bool getDragObjectActive() {
        return this->dragObjectActive;
    }

    void setDragObjectActive(bool set) {
        this->dragObjectActive = set;
    }

    bool getForceObjectActive() {
        return this->forceObjectActive;
    }

    void setForceObjectActive(bool active) {
        this->forceObjectActive = active;
    }

    bool getGeneratorActive() {
        return this->generatorActive;
    }

    void setGeneratorActive(bool active) {
        this->generatorActive = active;
    }

};