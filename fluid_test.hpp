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
#include "rendering.hpp"

class FluidHandler {
    int FLUID_CELL;
    int AIR_CELL;
    int SOLID_CELL;
    int numX;
    int numY;
    int n; // n = numY
    float invSpacing;
    float halfSpacing;
    float radius;
    float k;
    float dt;
    //std::vector<float> p;
    bool interacting = false;

    float checkseparationDist;
    float moveDist;

    float objectRadius;
    float objectX;
    float objectY;
    float objectPrevX;
    float objectPrevY;
    sf::CircleShape objectDrawer;
    float checkRigidObjectseparationDist;
    float objectXVel = 0;
    float objectYVel = 0;

    float mouseX = 0;
    float mouseY = 0;

    float forceObjectRadius = 250; // 200
    sf::CircleShape forceObjectDrawer;
    float checkForceObjectseparationDist;

    float generatorRadius = forceObjectRadius;
    sf::RectangleShape generatorDrawer;

    int numRowsPerThread;
    int numMissedRows;

    const float colorDiffusionCoeff = 0.001f;

    float diffusionRatio;

    float scalingFactor;
    int32_t scaledWIDTH;
    int32_t scaledHEIGHT;

    CollisionGrid collisionGrid;
    CollisionGrid cellOccupantsGrid;

    bool forceObjectActive = true;
    bool rigidObjectActive = false;
    bool generatorActive = false;

    float textureSizeX;
    float textureSizeY;

    float obstacleTextureSizeX;
    float obstacleTextureSizeY;

    float overRelaxation;
    float numPressureIters;

    std::vector<int32_t> nr0;
    std::vector<int32_t> nr1;
    std::vector<int32_t> nr2;
    std::vector<int32_t> nr3;

    std::vector<float> d0;
    std::vector<float> d1;
    std::vector<float> d2;
    std::vector<float> d3;

    int32_t renderPattern = 0;

    bool fireActive = false;

    std::vector<uint32_t> collisions;

    std::vector<double> residual;
    std::vector<double> pressure;
    std::vector<double> Adiag;
    std::vector<double> si;
    std::vector<double> li;
    std::vector<double> precon;
    std::vector<double> z;
    std::vector<double> search;

    sf::RectangleShape cellDrawer;

    uint32_t numThreads;
    uint32_t particlesPerThread;
    uint32_t numMissedParticles; 

    bool solidDrawing = false;

    sf::Font font;

    sf::Text text;

    bool leftMouseDown = false;
    bool rightMouseDown = false;

    std::vector<sf::Vector2i> obstaclePositions;
    sf::VertexArray obstacleVa{sf::PrimitiveType::Quads};

    sf::Texture obstacleTexture;

    sf::RenderStates obstacleStates;

    sf::RectangleShape pencil;

    int pencilRadius = 1;

    std::vector<double> dotProducts;

    float miscellaneousTime = 0.f;
    float CollisionTime = 0.f;
    float ObstacleCollisionTime = 0.f;
    float ToGridTime = 0.f;
    float DensityUpdateTime = 0.f;
    float ProjectionTime = 0.f;
    float ToParticlesTime = 0.f;
    float RenderingTime = 0.f;
    float FillGridTime = 0.f;
    float SimStepTime = 0.f;
    float steps = 0.f;
    float DotTime = 0.f;
    float EqualsPlusTime = 0.f;
    float preconditioningTime = 0.f;
    float matVecTime = 0.f;
    float scaledAddTime = 0.f;

    PressureSolver &pressure_solver;

    FluidState &fluid_attributes;

    FluidRenderer &fluid_renderer;

public:
    FluidHandler(float k, float overRelaxation_, float numPressureIters_, FluidState& fas, PressureSolver& ps, FluidRenderer& fr): k(k), overRelaxation(overRelaxation_), numPressureIters(numPressureIters_), fluid_attributes(fas), pressure_solver(ps), fluid_renderer(fr) {

            /*font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");
            text.setFont(font);
            text.setPosition(10, 10);
            text.setFillColor(sf::Color::White);*/

            
            FLUID_CELL = fluid_attributes.FLUID_CELL;
            AIR_CELL = fluid_attributes.AIR_CELL;
            SOLID_CELL = fluid_attributes.SOLID_CELL;

            numY = fluid_attributes.numY;
            numX = fluid_attributes.numX;
            n = numY;

            invSpacing =  1.f / fluid_attributes.cellSpacing;
            halfSpacing = 0.5f * fluid_attributes.cellSpacing;
            radius = fluid_attributes.radius;

            pencil.setSize(sf::Vector2f{fluid_attributes.cellSpacing, fluid_attributes.cellSpacing});
            pencil.setOrigin(halfSpacing, halfSpacing);
            pencil.setOutlineThickness(1);
            pencil.setOutlineColor(sf::Color::Black);

            numThreads = fluid_attributes.thread_pool.m_thread_count;
            particlesPerThread = fluid_attributes.num_particles / numThreads;
            numMissedParticles = fluid_attributes.num_particles - numThreads * particlesPerThread;

            dotProducts.resize(numThreads);

            this->Adiag.resize(numX * numY);
            this->si.resize(numX * numY);
            this->li.resize(numX * numY);
            this->precon.resize(numX * numY);
            this->z.resize(numX * numY);
            this->search.resize(numX * numY);
            this->residual.resize(numX * numY);
            this->pressure.resize(numX * numY);

            auto size = sf::Vector2f(fluid_attributes.cellSpacing, fluid_attributes.cellSpacing);
            this->cellDrawer.setSize(size);
            this->cellDrawer.setOutlineColor(sf::Color::White);
            this->cellDrawer.setOutlineThickness(1.f);

            this->collisions.resize(fluid_attributes.num_particles);

            this->nr0.resize(2 * fluid_attributes.num_particles);
            this->nr1.resize(2 * fluid_attributes.num_particles);
            this->nr2.resize(2 * fluid_attributes.num_particles);
            this->nr3.resize(2 * fluid_attributes.num_particles);

            this->d0.resize(2 * fluid_attributes.num_particles);
            this->d1.resize(2 * fluid_attributes.num_particles);
            this->d2.resize(2 * fluid_attributes.num_particles);
            this->d3.resize(2 * fluid_attributes.num_particles);

            size_t numObstacles = 2 * numX + 2 * (numY - 2);
            this->obstacleVa.resize(4 * numObstacles);
            this->obstaclePositions.resize(numObstacles);

            sf::Color gray = sf::Color(150, 150, 150);

            obstacleTexture.loadFromFile("gray_square.png");
            obstacleTexture.generateMipmap();
            auto const obstacle_texture_size = static_cast<sf::Vector2f>(obstacleTexture.getSize());
            for (int index = 0; index < numObstacles; ++index) {
                int i = 4 * index;
                obstacleVa[i].texCoords = {0.f, 0.f};
                obstacleVa[i + 1].texCoords = {obstacle_texture_size.x, 0.f};
                obstacleVa[i + 2].texCoords = {obstacle_texture_size.x, obstacle_texture_size.y};
                obstacleVa[i + 3].texCoords = {0.f, obstacle_texture_size.y};

                obstacleVa[i].color = gray;
                obstacleVa[i + 1].color = gray;
                obstacleVa[i + 2].color = gray;
                obstacleVa[i + 3].color = gray;
            }
            obstacleStates.texture = &obstacleTexture;

            obstacleTextureSizeX = obstacle_texture_size.x;
            obstacleTextureSizeY = obstacle_texture_size.y;

            this->scalingFactor = 2 * radius;

            this->scaledWIDTH = std::ceil(static_cast<float>(fluid_attributes.WIDTH) / scalingFactor);
            this->scaledHEIGHT = std::ceil(static_cast<float>(fluid_attributes.HEIGHT) / scalingFactor);

            collisionGrid = CollisionGrid(scaledWIDTH, scaledHEIGHT);

            cellOccupantsGrid = CollisionGrid(numX, numY);
            
            // move these to respective files
            int idx = 0;
            for (int i = 0; i < this->numX; ++i) {
                fluid_attributes.cellType[i * n] = SOLID_CELL;
                obstaclePositions[idx] = sf::Vector2i{i, 0};

                auto pos = gridCellToPos(i * n);
                float px = pos.x;
                float py = pos.y;
                size_t vaIdx = 4 * idx;
                obstacleVa[vaIdx].position = {px - halfSpacing, py - halfSpacing};
                obstacleVa[vaIdx + 1].position = {px + halfSpacing, py - halfSpacing};
                obstacleVa[vaIdx + 2].position = {px + halfSpacing, py + halfSpacing};
                obstacleVa[vaIdx + 3].position = {px - halfSpacing, py + halfSpacing};

                ++idx;
            }
            for (int i = 0; i < this->numX; ++i) {
                fluid_attributes.cellType[i * n + numY - 1] = SOLID_CELL;
                obstaclePositions[idx] = sf::Vector2i{i, numY - 1};

                auto pos = gridCellToPos(i * n + numY - 1);
                float px = pos.x;
                float py = pos.y;
                size_t vaIdx = 4 * idx;
                obstacleVa[vaIdx].position = {px - halfSpacing, py - halfSpacing};
                obstacleVa[vaIdx + 1].position = {px + halfSpacing, py - halfSpacing};
                obstacleVa[vaIdx + 2].position = {px + halfSpacing, py + halfSpacing};
                obstacleVa[vaIdx + 3].position = {px - halfSpacing, py + halfSpacing};

                ++idx;
            }
            for (int j = 1; j < this->numY - 1; ++j) {
                fluid_attributes.cellType[j] = SOLID_CELL;
                obstaclePositions[idx] = sf::Vector2i{0, j};

                auto pos = gridCellToPos(j);
                float px = pos.x;
                float py = pos.y;
                size_t vaIdx = 4 * idx;
                obstacleVa[vaIdx].position = {px - halfSpacing, py - halfSpacing};
                obstacleVa[vaIdx + 1].position = {px + halfSpacing, py - halfSpacing};
                obstacleVa[vaIdx + 2].position = {px + halfSpacing, py + halfSpacing};
                obstacleVa[vaIdx + 3].position = {px - halfSpacing, py + halfSpacing};

                ++idx;
            }
            for (int j = 1; j < this->numY - 1; ++j) {
                fluid_attributes.cellType[(numX - 1) * n + j] = SOLID_CELL;
                obstaclePositions[idx] = sf::Vector2i{numX - 1, j};

                auto pos = gridCellToPos((numX - 1) * n + j);
                float px = pos.x;
                float py = pos.y;
                size_t vaIdx = 4 * idx;
                obstacleVa[vaIdx].position = {px - halfSpacing, py - halfSpacing};
                obstacleVa[vaIdx + 1].position = {px + halfSpacing, py - halfSpacing};
                obstacleVa[vaIdx + 2].position = {px + halfSpacing, py + halfSpacing};
                obstacleVa[vaIdx + 3].position = {px - halfSpacing, py + halfSpacing};

                ++idx;
            }

            // -----------------------------------

            this->moveDist = 2 * radius;
            this->checkseparationDist = moveDist * moveDist;

            this->objectRadius = 50;//fluid_attributes.cellSpacing * 3;
            this->objectX = std::floor(fluid_attributes.WIDTH / 2);
            this->objectY = 2 * objectRadius + 10;
            this->objectPrevX = objectX;
            this->objectPrevY = objectY;
            this->objectDrawer.setOrigin(objectRadius, objectRadius);
            this->objectDrawer.setRadius(objectRadius);
            this->objectDrawer.setFillColor(sf::Color(255, 0, 0));
            this->checkRigidObjectseparationDist = (this->radius + objectRadius) * (this->radius + objectRadius);

            this->forceObjectDrawer.setOrigin(forceObjectRadius, forceObjectRadius);
            this->forceObjectDrawer.setRadius(forceObjectRadius);
            this->forceObjectDrawer.setOutlineThickness(1.f);
            this->forceObjectDrawer.setFillColor(sf::Color::Transparent);
            this->forceObjectDrawer.setOutlineColor(sf::Color::Red); 
            this->checkForceObjectseparationDist = (this->radius + forceObjectRadius) * (this->radius + forceObjectRadius);
            
            this->generatorDrawer.setOrigin(generatorRadius, generatorRadius);
            this->generatorDrawer.setSize(sf::Vector2f(2 * generatorRadius, 2 * generatorRadius));
            this->generatorDrawer.setOutlineThickness(1.f);
            this->generatorDrawer.setFillColor(sf::Color::Transparent);
            this->generatorDrawer.setOutlineColor(sf::Color::Red); 
    }

    sf::Vector2f gridCellToPos(int idx) {
        int i = idx % n;
        int j = idx / n;
        float x = halfSpacing + (j * fluid_attributes.cellSpacing);
        float y = halfSpacing + (i * fluid_attributes.cellSpacing);
        return sf::Vector2f{x, y};
    }

    // move all these into general_math file
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

    float clamp(float x, float min, float max) {
        return (x < min) * min + (x > max) * max + (x >= min && x <= max) * x;
    }

    float sign(float x) {
        return (x < 0) * -1.f + (x > 0) * 1.f;
    }

    void createRandomPositions() {
        std::uniform_int_distribution<int> randWidth(radius, fluid_attributes.WIDTH - radius);
        std::uniform_int_distribution<int> randHeight(radius, fluid_attributes.HEIGHT - radius);

        std::random_device rd;
        std::mt19937 mt(rd());

        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            fluid_attributes.positions[i * 2] = randWidth(mt);
            fluid_attributes.positions[i * 2 + 1] = randHeight(mt);
        }
    }

    void integrate(const uint32_t startIndex, const uint32_t endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            fluid_attributes.positions[2 * i] += fluid_attributes.velocities[2 * i] * dt;
            fluid_attributes.positions[2 * i + 1] += fluid_attributes.velocities[2 * i + 1] * dt;
            fluid_attributes.velocities[2 * i + 1] += fluid_attributes.gravityX * dt;
            fluid_attributes.velocities[2 * i] += fluid_attributes.gravityY * dt;

            if (this->fireActive && this->renderPattern == 3 && fluid_attributes.positions[2 * i + 1] < fluid_attributes.HEIGHT - fluid_attributes.cellSpacing - 10) {
                fluid_attributes.velocities[2 * i + 1] -= fluid_attributes.fireStrength * fluid_attributes.temperatures[i] * dt;
                if (fluid_attributes.temperatures[i] > 0) {
                    fluid_attributes.temperatures[i] -= fluid_attributes.tempDiffusion;
                }
                if (collisions[i] == 0) {
                    fluid_attributes.temperatures[i] = 0;
                }
            }
        }
    }

    void integrateMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->integrate(start, end);
        });
    }

    void addObjectsToGrids()
    {

        collisionGrid.clear();
        cellOccupantsGrid.clear();

        const float minX = fluid_attributes.cellSpacing;
        const float maxX = fluid_attributes.WIDTH - fluid_attributes.cellSpacing;
        const float minY = fluid_attributes.cellSpacing;
        const float maxY = fluid_attributes.HEIGHT - fluid_attributes.cellSpacing;

        uint32_t i{0};

        // makes it here then crashes
        for (int32_t index = 0; index < fluid_attributes.num_particles; ++index) {
            float x = fluid_attributes.positions[2 * index];
            float y = fluid_attributes.positions[2 * index + 1];
            if (x > minX && x < maxX && y > minY && y < maxY) {
                int32_t gridX = x / scalingFactor;
                int32_t gridY = y / scalingFactor;
                collisionGrid.addAtom(gridX, gridY, i);
                
                int32_t cellOccupantsX = x / fluid_attributes.cellSpacing;
                int32_t cellOccupantsY = y / fluid_attributes.cellSpacing;
                cellOccupantsGrid.addAtom(cellOccupantsX, cellOccupantsY, i);

            }
            ++i;
        }
    }

    void solveContact(uint32_t index, uint32_t otherIndex)
    {
        constexpr float eps           = 0.0001f;
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

            const float transfer = (fluid_attributes.temperatures[index] - fluid_attributes.temperatures[otherIndex]) * fluid_attributes.interConductivity * dt;
            fluid_attributes.temperatures[index] -= transfer;
            fluid_attributes.temperatures[otherIndex] += transfer;

            collisions[index]++;
            collisions[otherIndex]++;


        }
    }

    void checkAtomCellCollisions(uint32_t atom_idx, const CollisionCell& c)
    {
        for (uint32_t i{0}; i < c.objects_count; ++i) {
            solveContact(atom_idx, c.objects[i]);
        }
    }

    void processCell(const CollisionCell& c, uint32_t index)
    {
        for (uint32_t i{0}; i < c.objects_count; ++i) {
            const uint32_t atom_idx = c.objects[i];
            for (int32_t side = 0; side < 2; ++side) {
                checkAtomCellCollisions(atom_idx, collisionGrid.data[index + collisionGrid.height + side]);
            }
            for (int32_t side = 0; side < 2; ++side) {
                checkAtomCellCollisions(atom_idx, collisionGrid.data[index + side]);   
            }
            //checkAtomCellCollisions(atom_idx, grid.data[index - grid.height]);
            checkAtomCellCollisions(atom_idx, collisionGrid.data[index - collisionGrid.height + 1]);
        }
    }

    void solveCollisionThreaded(uint32_t start, uint32_t end)
    {
        for (uint32_t idx{start}; idx < end; ++idx) {
            processCell(collisionGrid.data[idx], idx);
        }
    }

    void makeForceObjectQueries(const int32_t strength) {
        const int32_t numCovered = std::ceil(forceObjectRadius / scalingFactor);

        const uint32_t mouseColumn = std::floor(mouseX / scalingFactor);
        const uint32_t mouseRow = std::floor(mouseY / scalingFactor);

        for (int32_t i = -numCovered; i < numCovered + 1; ++i) {
            for (int32_t j = -numCovered; j < numCovered + 1; ++j) {

                if (mouseRow + i <= 1 || mouseRow + i >= scaledHEIGHT - 1 || mouseColumn + j <= 1 || mouseColumn + j >= scaledWIDTH - 1) continue;

                const auto cell = collisionGrid.data[mouseRow + i + collisionGrid.height * (mouseColumn + j)];

                for (uint32_t i{0}; i < cell.objects_count; ++i) {
                    const uint32_t particleIndex = cell.objects[i];

                    float dx = fluid_attributes.positions[particleIndex * 2] - mouseX;
                    float dy = fluid_attributes.positions[particleIndex * 2 + 1] - mouseY;
                    float d2 = dx * dx + dy * dy;

                    if (d2 > checkForceObjectseparationDist || d2 == 0.0f) continue;

                    float d = std::sqrt(d2);

                    float edgeT = d / forceObjectRadius;
            
                    float centerT = 1 - edgeT;

                    fluid_attributes.velocities[2 * particleIndex] += (dx * strength - fluid_attributes.velocities[2 * particleIndex]) * centerT * dt;
                    fluid_attributes.velocities[2 * particleIndex + 1] += (dy * strength - fluid_attributes.velocities[2 * particleIndex + 1]) * centerT * dt;
                }
            }
        }
    }

    void solveCollisions()
    {
        const uint32_t slice_count  = numThreads * 2;
        const uint32_t slice_size   = (collisionGrid.width / slice_count) * collisionGrid.height;
        const uint32_t last_cell    = 2 * numThreads * slice_size;
        
        for (uint32_t i{0}; i < numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([this, i, slice_size]{
                uint32_t const start{2 * i * slice_size};
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
        
        for (uint32_t i{0}; i < numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([this, i, slice_size]{
                uint32_t const start{(2 * i + 1) * slice_size};
                uint32_t const end  {start + slice_size};
                solveCollisionThreaded(start, end);
            });
        }
        fluid_attributes.thread_pool.waitForCompletion();
    }

    void constrainWalls(const uint32_t startIndex, const uint32_t endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            if (fluid_attributes.positions[2 * i] - radius < this->fluid_attributes.cellSpacing) {
                fluid_attributes.positions[2 * i] = radius + this->fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[2 * i] < 0) {
                    fluid_attributes.velocities[2 * i] = 0.f;
                }
            }
            else if (fluid_attributes.positions[2 * i] + radius > fluid_attributes.WIDTH - this->fluid_attributes.cellSpacing) {
                fluid_attributes.positions[2 * i] = fluid_attributes.WIDTH - radius - this->fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[2 * i] > 0) {
                    fluid_attributes.velocities[2 * i] = 0.f;
                }
            }
            if (fluid_attributes.positions[2 * i + 1] - radius < this->fluid_attributes.cellSpacing) {
                fluid_attributes.positions[2 * i + 1] = radius + this->fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[2 * i + 1] < 0) {
                    fluid_attributes.velocities[2 * i + 1] = 0.f;
                }
            }
            else if (fluid_attributes.positions[2 * i + 1] + radius > fluid_attributes.HEIGHT - this->fluid_attributes.cellSpacing) {
                fluid_attributes.positions[2 * i + 1] = fluid_attributes.HEIGHT - radius - this->fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[2 * i + 1] > 0) {
                    fluid_attributes.velocities[2 * i + 1] = 0.f;
                }
                const float remove = 10.f; // 5
                if (this->renderPattern == 3 && fluid_attributes.temperatures[i] < fluid_renderer.tempgradient.size() && fluid_attributes.positions[2 * i] > fluid_attributes.WIDTH / remove && fluid_attributes.positions[2 * i] < fluid_attributes.WIDTH - fluid_attributes.WIDTH / remove) {
                    fluid_attributes.temperatures[i] += fluid_attributes.groundConductivity;
                }
            }
        }
    }

    void constrainWallsMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->constrainWalls(start, end);
        });
    }

    float calculateBoxNormals(float bx, float by, float px, float py, float& nx, float& ny) {
        float localpx = px - bx;
        float localpy = py - by;
        float dx = std::abs(localpx) - halfSpacing;
        float dy = std::abs(localpy) - halfSpacing;
        float pdx = std::max(dx, 0.f);
        float pdy = std::max(dy, 0.f);
        float insideDist = std::min(std::max(dx, dy), 0.f);
        float dist = sqrt(pdx * pdx + pdy * pdy) + insideDist - radius;
    
        if (dist >= 0.f) {
            return dist;
        }
    
        float dirX = sign(localpx);
        float dirY = sign(localpy);

        int idx = n * static_cast<int>(px * invSpacing) + static_cast<int>(py * invSpacing);
    
        if (pdx > 0 && pdy > 0 && fluid_attributes.cellType[idx + dirX] != SOLID_CELL && fluid_attributes.cellType[idx + dirY] != SOLID_CELL) {
            float len = dist - (insideDist - radius);
            nx = -dirX * dx / len;
            ny = -dirY * dy / len;
        }
        else if (abs(localpx) > abs(localpy)) {
            nx = dirX * sign(dist);
            ny = 0;
        }
        else {
            nx = 0;
            ny = dirY * sign(dist);
        }
    
        return dist;
    }

    void separateParticle(float& px, float& py, float& vx, float& vy) {
        int localX = px / fluid_attributes.cellSpacing;
        int localY = py / fluid_attributes.cellSpacing;

        int x0 = std::max(0, localX - 1);
        int x1 = std::min(numX - 1, localX + 1);
        int y0 = std::max(0, localY - 1);
        int y1 = std::min(numY - 1, localY + 1);

        float prevX = px;
        float prevY = py;

        for (int i = x0; i <= x1; ++i) {
            for (int j = y0; j <= y1; ++j) {
                if (fluid_attributes.cellType[i * n + j] == SOLID_CELL) {
                    float nx = 0;
                    float ny = 0;
                    float dist = calculateBoxNormals(i * fluid_attributes.cellSpacing + halfSpacing, j * fluid_attributes.cellSpacing + halfSpacing, px, py, nx, ny);
                    if (dist < 0.f) {
                        float velocityNormal = vx * -nx + vy * -ny;
                        if (velocityNormal < 0) {
                            vx -= velocityNormal * -nx;
                            vy -= velocityNormal * -ny;
                        }
                        px += nx * dist;
                        py += ny * dist;
                    }
                }
            }
        }
    }

    void collideSurfaces(const uint32_t start, const uint32_t end) {
        for (int i = start; i < end; ++i) {
            separateParticle(fluid_attributes.positions[2 * i], fluid_attributes.positions[2 * i + 1], fluid_attributes.velocities[2 * i], fluid_attributes.velocities[2 * i + 1]);
        }
    }

    void collideSurfacesMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->collideSurfaces(start, end);
        });
    }

    void drawSolids() {
        int localX = static_cast<int>(mouseX / fluid_attributes.cellSpacing);
        int localY = static_cast<int>(mouseY / fluid_attributes.cellSpacing);

        int numObstacles = obstaclePositions.size();

        int numPotentialAddedObstacles = ((2 * pencilRadius + 1) * (2 * pencilRadius + 1));
        obstaclePositions.resize(numObstacles + numPotentialAddedObstacles);

        int numAddedObstacles = 0;
        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                int x = localX + i;
                int y = localY + j;

                if (x % numX > 0 && y > 0 && x % numX < numX - 1 && y < numY - 1) {

                    int idx = x * numY + y;

                    if (fluid_attributes.cellType[idx] != SOLID_CELL) {
                        fluid_attributes.cellType[idx] = SOLID_CELL;
                        obstaclePositions[numObstacles + numAddedObstacles] = sf::Vector2i{x, y};
                        numAddedObstacles++;
                    }
                }
            }
        }

        obstaclePositions.resize(numObstacles + numAddedObstacles);

        if (numAddedObstacles > 0) {
            obstacleVa.resize(4 * (numObstacles + numAddedObstacles));

            sf::Color gray = sf::Color(150, 150, 150);
            for (int i = 0; i < numAddedObstacles; ++i) {
                int idx = 4 * (i + numObstacles);
                obstacleVa[idx].texCoords = {0.f, 0.f};
                obstacleVa[idx + 1].texCoords = {textureSizeX, 0.f};
                obstacleVa[idx + 2].texCoords = {textureSizeX, textureSizeY};
                obstacleVa[idx + 3].texCoords = {0.f, textureSizeY};

                obstacleVa[idx].color = gray;
                obstacleVa[idx + 1].color = gray;
                obstacleVa[idx + 2].color = gray;
                obstacleVa[idx + 3].color = gray;

                auto cellCoords = obstaclePositions[numObstacles + i];
                auto pos = gridCellToPos(cellCoords.x * n + cellCoords.y);
                float px = pos.x;
                float py = pos.y;

                obstacleVa[idx].position = {px - halfSpacing, py - halfSpacing};
                obstacleVa[idx + 1].position = {px + halfSpacing, py - halfSpacing};
                obstacleVa[idx + 2].position = {px + halfSpacing, py + halfSpacing};
                obstacleVa[idx + 3].position = {px - halfSpacing, py + halfSpacing};
            }
        }    
    }

    void eraseSolids() {
        int localX = static_cast<int>(mouseX / fluid_attributes.cellSpacing);
        int localY = static_cast<int>(mouseY / fluid_attributes.cellSpacing);

        int numObstacles = obstaclePositions.size();
        int numFreedCells = 0;
        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                int x = localX + i;
                int y = localY + j;

                int idx = x * n + y;
                if (fluid_attributes.cellType[idx] == SOLID_CELL && x > 0 && y > 0 && x < numX - 1 && y < numY - 1) {
                    fluid_attributes.cellType[idx] = AIR_CELL;
                    // find the index of obstacleSet with the position we want to remove
                    // replace that index with the [numObstacles - (numFreedCells + 1)] index of obstacleSet
                    // then just resize the obstacleSet and obstacleVa after this double loop
                    numFreedCells++;

                    int obstacleIdx = find(&obstaclePositions, sf::Vector2i{x, y});
                    int vaIdx = 4 * obstacleIdx;

                    int posMv = numObstacles - numFreedCells;
                    int vaMv = 4 * posMv;

                    obstaclePositions[obstacleIdx] = obstaclePositions[posMv];

                    obstacleVa[vaIdx].position = obstacleVa[vaMv].position;
                    obstacleVa[vaIdx + 1].position = obstacleVa[vaMv + 1].position;
                    obstacleVa[vaIdx + 2].position = obstacleVa[vaMv + 2].position;
                    obstacleVa[vaIdx + 3].position = obstacleVa[vaMv + 3].position;
                }
            }
        }
        obstacleVa.resize(4 * (numObstacles - numFreedCells));
        obstaclePositions.resize(numObstacles - numFreedCells);
    }

    void cacheTransferNodes(int32_t start, int32_t end, float halfHeight, int32_t component) {
        const float h2 = halfHeight;

        const float dx = (component != 0) * h2;
        const float dy = (component == 0) * h2;

        for (int32_t i = start; i < end; ++i) {
            float x = fluid_attributes.positions[2 * i];
            float y = fluid_attributes.positions[2 * i + 1];
            x = this->clamp(x, fluid_attributes.cellSpacing, (this->numX - 1) * fluid_attributes.cellSpacing);
            y = this->clamp(y, fluid_attributes.cellSpacing, (this->numY - 1) * fluid_attributes.cellSpacing);
            // x0 is the grid position to the left of the particle, x1 is the position to the right of the particle. Both can only go up to the second to last cell to the right in the grid because we dont want to be changing wall velocities
            int x0 = std::max(1, std::min(static_cast<int>(std::floor((x - dx) * invSpacing)), this->numX - 2)); // - 1
            // basically x - xCell to get the weight of that cell in relation to the particle
            // in this case, x is moved over to x - dx, and xCell is just grid position of x multiplied by grid spacing
            float tx = ((x - dx) - x0 * fluid_attributes.cellSpacing) * invSpacing;
            // add 1 to get the cell to the right
            int x1 = std::min(x0 + 1, this->numX - 2); // - 1
            // this fixes a bug that makes water touching the left wall and ceiling explode sometimes 
            if (component == 0 && x0 == 1) {
                x0 = x1;
            }
            if (component == 1 && x0 == 1) {
                x1 = x0;
            }
            // same thing with y
            int y0 = std::max(0, std::min(static_cast<int>(std::floor((y - dy) * invSpacing)), this->numY - 2)); // - 2
            float ty = ((y - dy) - y0 * fluid_attributes.cellSpacing) * invSpacing;
            int y1 = std::min(y0 + 1, this->numY - 1); // - 1
            // try seeing if these have any problems when CG is implemented
            // fixes jitter when the fluid is held with the force object at the top, but also makes the fluid stick to the top so leave out. Could just only do this first one when gravity < 0
            /*if (component == 1 && y0 == 0) {
                y0 = y1;
            }
            if (component == 0 && y0 == 0) {
                y1 = y0;
            }*/
            // fixes jitter when the fluid is held with the force object at the bottom, but also introduces the same gauss seidel artifact as the right side
            /*if (component == 1 && y0 == numY - 2) {
                y1 = y0;
            }*/
            float sx = 1.f - tx;
            float sy = 1.f - ty;
            // weights for each corner in u/v field
            d0[2 * i + component] = sx * sy;
            d1[2 * i + component] = tx * sy;
            d2[2 * i + component] = tx * ty;
            d3[2 * i + component] = sx * ty;

            // top left
            nr0[2 * i + component] = x0 * n + y0;
            // top right
            nr1[2 * i + component] = x1 * n + y0;
            // bottom right
            nr2[2 * i + component] = x1 * n + y1;
            //bottom left
            nr3[2 * i + component] = x0 * n + y1;
        }
    }

    void cacheTransferNodesMulti() {
        for (int32_t i = 0; i < numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([&, i](){
                cacheTransferNodes(fluid_attributes.particlesPerThread * i, fluid_attributes.particlesPerThread * i + fluid_attributes.particlesPerThread, halfSpacing, 0);
                cacheTransferNodes(fluid_attributes.particlesPerThread * i, fluid_attributes.particlesPerThread * i + fluid_attributes.particlesPerThread, halfSpacing, 1);
            });
        }

        cacheTransferNodes(fluid_attributes.num_particles - numMissedParticles, fluid_attributes.num_particles, halfSpacing, 0);
        cacheTransferNodes(fluid_attributes.num_particles - numMissedParticles, fluid_attributes.num_particles, halfSpacing, 1);

        fluid_attributes.thread_pool.waitForCompletion();
    }

    void transfer(const bool& toGrid) {
        if (toGrid) {
            this->cacheTransferNodesMulti();
            this->setUpTransferGrids();
            this->transferToGrid();
        }
        else {
            this->transferToParticles();
        }
    }

    void setUpTransferGrids() {
        const float n = this->numY;
        const float h2 = 0.5 * fluid_attributes.cellSpacing;

        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));
        std::fill(begin(fluid_attributes.du), end(fluid_attributes.du), 0.f);
        std::fill(begin(fluid_attributes.dv), end(fluid_attributes.dv), 0.f);
        std::fill(begin(fluid_attributes.u), end(fluid_attributes.u), 0.f);
        std::fill(begin(fluid_attributes.v), end(fluid_attributes.v), 0.f);

        // initialize every inside cell to air
        for (int i = 0; i < this->numX; ++i) {
            for (int j = 0; j < this->numY; ++j) {
                if (fluid_attributes.cellType[i * n + j] != SOLID_CELL) {
                    fluid_attributes.cellType[i * n + j] = AIR_CELL;
                }
            }
        }

        // initialize all cells that particles are in to fluid
        for (int i = 0; i < this->fluid_attributes.num_particles; ++i) {
            float x = fluid_attributes.positions[2 * i];
            float y = fluid_attributes.positions[2 * i + 1];

            int xi = this->clamp(std::floor(x * invSpacing), 0, this->numX - 1);
            int yi = this->clamp(std::floor(y * invSpacing), 0, this->numY - 1);

            int cellNr = xi * n + yi;
            if (fluid_attributes.cellType[cellNr] == AIR_CELL) {
                fluid_attributes.cellType[cellNr] = FLUID_CELL;
            }
        }
    }

    void transferToParticle(int32_t startIndex, int32_t endIndex) {
        const int32_t n = this->numY;
        for (int32_t i = startIndex; i < endIndex; ++i) {
            const int32_t ui = 2 * i;
            const int32_t vi = 2 * i + 1;
            const int32_t nr0_u = nr0[ui];
            const int32_t nr1_u = nr1[ui];
            const int32_t nr2_u = nr2[ui];
            const int32_t nr3_u = nr3[ui];

            const float d0_u = d0[ui];
            const float d1_u = d1[ui];
            const float d2_u = d2[ui];
            const float d3_u = d3[ui];

            const int32_t nr0_v = nr0[vi];
            const int32_t nr1_v = nr1[vi];
            const int32_t nr2_v = nr2[vi];
            const int32_t nr3_v = nr3[vi];

            const float d0_v = d0[vi];
            const float d1_v = d1[vi];
            const float d2_v = d2[vi];
            const float d3_v = d3[vi];

            const float pvx = fluid_attributes.velocities[ui];
            const float pvy = fluid_attributes.velocities[vi];
           
            // these will be used to make sure that air cells are not considered when transferring velocities back to particles
            // nr0 - offset is the same thing as [(i-1) * n]
            // if u is being considered, then we only have to check left and right cells ([nr0] and [nr0 - n])
            // if v is being considered, then we only have to check above and below cells ([nr0] and [nr0 - 1])
            const float valid0u = fluid_attributes.cellType[nr0_u] != AIR_CELL || fluid_attributes.cellType[nr0_u - n] != AIR_CELL;
            const float valid1u = fluid_attributes.cellType[nr1_u] != AIR_CELL || fluid_attributes.cellType[nr1_u - n] != AIR_CELL;
            const float valid2u = fluid_attributes.cellType[nr2_u] != AIR_CELL || fluid_attributes.cellType[nr2_u - n] != AIR_CELL;
            const float valid3u = fluid_attributes.cellType[nr3_u] != AIR_CELL || fluid_attributes.cellType[nr3_u - n] != AIR_CELL;

            const float valid0v = fluid_attributes.cellType[nr0_v] != AIR_CELL || fluid_attributes.cellType[nr0_v - 1] != AIR_CELL;
            const float valid1v = fluid_attributes.cellType[nr1_v] != AIR_CELL || fluid_attributes.cellType[nr1_v - 1] != AIR_CELL;
            const float valid2v = fluid_attributes.cellType[nr2_v] != AIR_CELL || fluid_attributes.cellType[nr2_v - 1] != AIR_CELL;
            const float valid3v = fluid_attributes.cellType[nr3_v] != AIR_CELL || fluid_attributes.cellType[nr3_v - 1] != AIR_CELL;

            const float divX = valid0u * d0_u + valid1u * d1_u + valid2u * d2_u + valid3u * d3_u;
            const float divY = valid0v * d0_v + valid1v * d1_v + valid2v * d2_v + valid3v * d3_v;

            float picV;
            float corr;
            float flipV;

            if (divX > 0.f) {
                picV = (valid0u * d0_u * fluid_attributes.u[nr0_u] + valid1u * d1_u * fluid_attributes.u[nr1_u] + valid2u * d2_u * fluid_attributes.u[nr2_u] + valid3u * d3_u * fluid_attributes.u[nr3_u]) / divX;

                corr = (valid0u * d0_u * (fluid_attributes.u[nr0_u] - fluid_attributes.prevU[nr0_u]) + valid1u * d1_u * (fluid_attributes.u[nr1_u] - fluid_attributes.prevU[nr1_u]) + valid2u * d2_u * (fluid_attributes.u[nr2_u] - fluid_attributes.prevU[nr2_u]) + valid3u * d3_u * (fluid_attributes.u[nr3_u] - fluid_attributes.prevU[nr3_u])) / divX;
                flipV = pvx + corr;
                fluid_attributes.velocities[ui] = (1.f - fluid_attributes.flipRatio) * picV + fluid_attributes.flipRatio * flipV;
            }

            if (divY > 0.f) {
                picV = (valid0v * d0_v * fluid_attributes.v[nr0_v] + valid1v * d1_v * fluid_attributes.v[nr1_v] + valid2v * d2_v * fluid_attributes.v[nr2_v] + valid3v * d3_v * fluid_attributes.v[nr3_v]) / divY;

                corr = (valid0v * d0_v * (fluid_attributes.v[nr0_v] - fluid_attributes.prevV[nr0_v]) + valid1v * d1_v * (fluid_attributes.v[nr1_v] - fluid_attributes.prevV[nr1_v]) + valid2v * d2_v * (fluid_attributes.v[nr2_v] - fluid_attributes.prevV[nr2_v]) + valid3v * d3_v * (fluid_attributes.v[nr3_v] - fluid_attributes.prevV[nr3_v])) / divY;
                flipV = pvy + corr;
                fluid_attributes.velocities[vi] = (1.f - fluid_attributes.flipRatio) * picV + fluid_attributes.flipRatio * flipV;
            }
        }
    }

    void transferToParticles() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->transferToParticle(start, end);
        });
    }

    void transferToUGridCells(int idx) {
        const auto cell = cellOccupantsGrid.data[idx];
    
        for (uint32_t id{0}; id < cell.objects_count; ++id) {
            const uint32_t particleIndex = cell.objects[id];

            const int32_t ui = 2 * particleIndex;
            const int32_t vi = 2 * particleIndex + 1;

            const int32_t nr0_u = nr0[ui];
            const int32_t nr1_u = nr1[ui];
            const int32_t nr2_u = nr2[ui];
            const int32_t nr3_u = nr3[ui];

            const float d0_u = d0[ui];
            const float d1_u = d1[ui];
            const float d2_u = d2[ui];
            const float d3_u = d3[ui];

            const float pvx = fluid_attributes.velocities[ui];
            const float pvy = fluid_attributes.velocities[vi];
           
            fluid_attributes.u[nr0_u] += pvx * d0_u;  
            fluid_attributes.u[nr1_u] += pvx * d1_u;
            fluid_attributes.u[nr2_u] += pvx * d2_u;
            fluid_attributes.u[nr3_u] += pvx * d3_u;

            fluid_attributes.du[nr0_u] += d0_u;
            fluid_attributes.du[nr1_u] += d1_u;
            fluid_attributes.du[nr2_u] += d2_u;
            fluid_attributes.du[nr3_u] += d3_u;
        }
    }

    void transferToUGrid(int start, int end) {

        for (int i = start; i < end; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int idx = i * n + j;
                transferToUGridCells(idx);
            }
        }
    }

    void transferToVGridCells(int idx) {
        const auto cell = cellOccupantsGrid.data[idx];
    
        for (uint32_t id{0}; id < cell.objects_count; ++id) {
            const uint32_t particleIndex = cell.objects[id];

            const int32_t ui = 2 * particleIndex;
            const int32_t vi = 2 * particleIndex + 1;

            const int32_t nr0_v = nr0[vi];
            const int32_t nr1_v = nr1[vi];
            const int32_t nr2_v = nr2[vi];
            const int32_t nr3_v = nr3[vi];

            const float d0_v = d0[vi];
            const float d1_v = d1[vi];
            const float d2_v = d2[vi];
            const float d3_v = d3[vi];

            const float pvx = fluid_attributes.velocities[ui];
            const float pvy = fluid_attributes.velocities[vi];
    
            fluid_attributes.v[nr0_v] += pvy * d0_v;  
            fluid_attributes.v[nr1_v] += pvy * d1_v;
            fluid_attributes.v[nr2_v] += pvy * d2_v;
            fluid_attributes.v[nr3_v] += pvy * d3_v;

            fluid_attributes.dv[nr0_v] += d0_v;
            fluid_attributes.dv[nr1_v] += d1_v;
            fluid_attributes.dv[nr2_v] += d2_v;
            fluid_attributes.dv[nr3_v] += d3_v;
        }
    }

    void transferToVGrid(int start, int end) {

        for (int i = start; i < end; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int idx = i * n + j;
                transferToVGridCells(idx);
            }
        }
    }

    void transferToGrid() {

        const int32_t halfNumRemainingThreads = numThreads > 3 ? (numThreads - 2) / 2 : 1;
        const int32_t numColumnsPerThread = (numX - 2) / halfNumRemainingThreads;
        const int32_t numMissedColumns = numX - 2 - numColumnsPerThread * halfNumRemainingThreads;

        for (int i = 0; i < halfNumRemainingThreads; ++i) {
            int start = i * numColumnsPerThread;
            int end = (i == halfNumRemainingThreads - 1) ? (numX - 1) : (start + numColumnsPerThread);
            fluid_attributes.thread_pool.addTask([&, start, end]() {
                this->transferToUGrid(start, end);
            });
            fluid_attributes.thread_pool.addTask([&, start, end]() {
                this->transferToVGrid(start, end);
            });
        }

        fluid_attributes.thread_pool.waitForCompletion();


        for (int i = 0; i < fluid_attributes.u.size(); ++i) {
            float prevNode = fluid_attributes.du[i];
            if (prevNode > 0.f) {
                fluid_attributes.u[i] /= prevNode;
            }
        }
        for (int i = 0; i < fluid_attributes.v.size(); ++i) {
            float prevNode = fluid_attributes.dv[i];
            if (prevNode > 0.f) {
                fluid_attributes.v[i] /= prevNode;
            }
        }

        for (int i = 0; i < numX; ++i) {
            for (int j = 0; j < numY; ++j) {
                int idx = i * n + j;
                bool solid = fluid_attributes.cellType[idx] == SOLID_CELL;
                if (solid || i > 0 && fluid_attributes.cellType[idx - n] == SOLID_CELL) {
                    fluid_attributes.u[idx] = fluid_attributes.prevU[idx];
                }
                if (solid || j > 0 && fluid_attributes.cellType[idx - 1] == SOLID_CELL) {
                    fluid_attributes.v[idx] = fluid_attributes.prevV[idx];
                }
            }
        }
    }

    void updateParticleDensity() {

        std::fill(begin(fluid_attributes.cellDensities), end(fluid_attributes.cellDensities), 0.f);

        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            float x = fluid_attributes.positions[2 * i];
            float y = fluid_attributes.positions[2 * i + 1];

            x = this->clamp(x, fluid_attributes.cellSpacing, (this->numX - 1) * fluid_attributes.cellSpacing);
            y = this->clamp(y, fluid_attributes.cellSpacing, (this->numY - 1) * fluid_attributes.cellSpacing);

            int x0 = std::max(1, std::min((int)(std::floor((x - halfSpacing) * invSpacing)), this->numX - 2));
            float tx = ((x - halfSpacing) - x0 * fluid_attributes.cellSpacing) * invSpacing;
            int x1 = std::min(x0 + 1, this->numX - 1);

            int y0 = std::max(1, std::min((int)(std::floor((y - halfSpacing) * invSpacing)), this->numY - 2));
            float ty = ((y - halfSpacing) - y0 * fluid_attributes.cellSpacing) * invSpacing;
            int y1 = std::min(y0 + 1, this->numY - 2);
           
            float sx = 1.f - tx;
            float sy = 1.f - ty;

            if (x0 < this->numX && y0 < this->numY) {
                fluid_attributes.cellDensities[x0 * n + y0] += sx * sy;
            }
            if (x1 < this->numX && y0 < this->numY) {
                fluid_attributes.cellDensities[x1 * n + y0] += tx * sy;
            }
            if (x1 < this->numX && y1 < this->numY) {
                fluid_attributes.cellDensities[x1 * n + y1] += tx * ty;
            }
            if (x0 < this->numX && y1 < this->numY) {
                fluid_attributes.cellDensities[x0 * n + y1] += sx * ty;
            }
        }

        if (fluid_attributes.particleRestDensity == 0.f) {
            float sum = 0.f;
            int numFluidCells = 0;

            for (int i = 0; i < fluid_attributes.gridSize; ++i) {
                if (fluid_attributes.cellType[i] == FLUID_CELL) {
                    sum += fluid_attributes.cellDensities[i];
                    numFluidCells++;
                }
            }

            if (numFluidCells > 0) {
                fluid_attributes.particleRestDensity = sum / numFluidCells;
            }
        }
    }

    void setUpResidual() {
        double scale = 1.0 / fluid_attributes.cellSpacing;

        std::fill(begin(residual), end(residual), 0.0);

        for (int32_t i = 1; i < numX - 1; ++i) {
            for (int32_t j = 1; j < numY - 1; ++j) {
                int32_t idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID_CELL)  {
                    continue;
                }

                float divergence = -(fluid_attributes.u[idx + n] - fluid_attributes.u[idx] + fluid_attributes.v[idx + 1] - fluid_attributes.v[idx]);

                if (fluid_attributes.particleRestDensity > 0.f) {
                    float compression = fluid_attributes.cellDensities[idx] - fluid_attributes.particleRestDensity;
                    divergence += (compression > 0.f) * this->k * compression;
                }

                residual[idx] = divergence;

            }
        }
    }

    void ScaledAdd(std::vector<double>& a, std::vector<double>& b, double c) {
        for (int i = 0; i < numX * numY; ++i) {
            if (fluid_attributes.cellType[i] == FLUID_CELL) {
                a[i] += b[i] * c;
            }
        }
    }

    double DotMulti(std::vector<double> *a, std::vector<double> *b) {
        const int32_t numThreads_ = 1;
        const int32_t numColumnsPerThread = (numX * numY) / numThreads_;
        const int32_t numMissedColumns = numX * numY - numColumnsPerThread * numThreads_;

        std::fill(begin(dotProducts), end(dotProducts), 0.0);

        for (int i = 0; i < numThreads_; ++i) {
            int start = i * numColumnsPerThread;
            int end = (i == numThreads_ - 1) ? (numX * numY) : (start + numColumnsPerThread);
            fluid_attributes.thread_pool.addTask([&, i, start, end]() {
                this->Dot(a, b, start, end, dotProducts[i]);
            });
        }

        fluid_attributes.thread_pool.waitForCompletion();

        double res = 0.0;
        for (double el : dotProducts) {
            res += el;
        }

        return res;
    }

    void Dot(std::vector<double> *a, std::vector<double> *b, int start, int end, double& res) {
        for (int i = start; i < end; ++i) {
            if (fluid_attributes.cellType[i] == FLUID_CELL) {
                res += (*a)[i] * (*b)[i];
            }
        }
    }

    void EqualsPlusTimesMulti(std::vector<double> *a, std::vector<double> *b, double c) {
        const int32_t numThreads_ = 1;
        const int32_t numColumnsPerThread = (numX * numY) / numThreads_;
        const int32_t numMissedColumns = numX * numY - numColumnsPerThread * numThreads_;

        for (int i = 0; i < numThreads_; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->EqualsPlusTimes(a, b, c, i * numColumnsPerThread, i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->EqualsPlusTimes(a, b, c, numX * numY - numMissedColumns, numX * numY);

        fluid_attributes.thread_pool.waitForCompletion();

        /*const int32_t numColumnsPerThread = (numX * numY) / numThreads_;
        fluid_attributes.thread_pool.dispatch(numColumns, [this](int start, int end) {
            this->updateVertexArrayVorticity(start, end);
        });*/
    }

    void EqualsPlusTimes(std::vector<double> *a, std::vector<double> *b, double c, int start, int end) {
        for (int i = start; i < end; ++i) {
            if (fluid_attributes.cellType[i] == FLUID_CELL) {
                (*a)[i] = (*b)[i] + (*a)[i] * c;
            }
        }
    }

    void applyPressureMulti() {
        const int32_t numThreads_ = 1;
        const int32_t numColumnsPerThread = (numX - 2) / numThreads_;
        const int32_t numMissedColumns = (numX - 2) - numColumnsPerThread * numThreads_;

        for (int i = 0; i < numThreads_; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->applyPressure(1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->applyPressure(numX - 1 - numMissedColumns, numX - 1);

        fluid_attributes.thread_pool.waitForCompletion();
    }

    void applyPressure(int start, int end) {
        //const float density = 1000.f;
        //const float scale = dt / (density * fluid_attributes.cellSpacing);
    
        for (int i = start; i < end; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int idx = i * n + j;
                int leftIdx =  idx - n;
                int upIdx = idx - 1;

                if ((fluid_attributes.cellType[idx] == FLUID_CELL || fluid_attributes.cellType[leftIdx] == FLUID_CELL) &&
                    fluid_attributes.cellType[idx] != SOLID_CELL && fluid_attributes.cellType[leftIdx] != SOLID_CELL) {
                    float p = pressure[idx];
                    float pLeft = pressure[leftIdx];
                    fluid_attributes.u[idx] -= 1 * (p - pLeft);
                }
    
                if ((fluid_attributes.cellType[idx] == FLUID_CELL || fluid_attributes.cellType[upIdx] == FLUID_CELL) &&
                    fluid_attributes.cellType[idx] != SOLID_CELL && fluid_attributes.cellType[upIdx] != SOLID_CELL) {
                    float p = pressure[idx];
                    float pTop = pressure[upIdx];
                    fluid_attributes.v[idx] -= 1 * (p - pTop);
                }
            }
        }
    }

    void setUpA() {
        float scale = 1;//dt / (fluid_attributes.cellSpacing * fluid_attributes.cellSpacing); // also see what happens when you divide by 1000 (density of water) as well

        std::fill(Adiag.begin(), Adiag.end(), 0.0);

        // Ax[i] is the cell to the right, Ay[i] is the cell below
        std::fill(si.begin(), si.end(), 0.0);
        std::fill(li.begin(), li.end(), 0.0);

        for (int i = 1; i < numX - 1; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int idx = i * n + j;

                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;

                int left = fluid_attributes.cellType[idx - n];
                int right = fluid_attributes.cellType[idx + n];
                int bottom = fluid_attributes.cellType[idx + 1];
                int top = fluid_attributes.cellType[idx - 1];

                if (left == FLUID_CELL || left == AIR_CELL) { 
                    Adiag[idx] += scale;
                }

                if (right == FLUID_CELL) {
                    Adiag[idx] += scale;
                    li[idx] = -scale;
                }
                else if (right == AIR_CELL) {
                    Adiag[idx] += scale;
                }

                if (top == FLUID_CELL || top == AIR_CELL) { 
                    Adiag[idx] += scale;
                }

                if (bottom == FLUID_CELL) {
                    Adiag[idx] += scale;
                    si[idx] = -scale;
                }
                else if (bottom == AIR_CELL) {
                    Adiag[idx] += scale;
                }
            }
        }
    }

    void buildPreconditioner() {
        const double tau = 0.97;
        const double sigma = 0.25;

        for (int i = 1; i < numX - 1; ++i) { // I call this the level iteration. Every time j increases, we are at the next level iteration
            for (int j = 1; j < numY - 1; ++j) { // and this is the swipe iteration. Every time i increases, we are at the next swipe iteration
                int idx = i * n + j;

                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;

                double e = Adiag[idx];

                // si is short for swipe iteration, just some slang for you new gens
                    // si[index] = 0 if the cell in the next swipe iteration touching the current is not fluid, else si[index] = -1

                // li is short for level iteration, try and keep up with the slang
                    // li[index] = 0 if the cell in the next level iteration touching the current cell is not fluid, else li[index] = -1

                // if I add/subtract something by S in a comment, that means I'm taking the value of it at the next/previous swipe iteration
                // same for level iteration, denoted with L

                // if cell - S is fluid, then do: (si - S) * (precon - S), and (li - S) * (precon - S)
                if (fluid_attributes.cellType[idx - 1] == FLUID_CELL) {
                    double px = si[idx - 1] * precon[idx - 1];
                    double py = li[idx - 1] * precon[idx - 1];
                    e -= (px * px + tau * px * py);
                }

                // if cell - L is fluid, then do: (si - L) * (precon - L), and (li - L) * (precon - L)
                if (fluid_attributes.cellType[idx - n] == FLUID_CELL) {
                    double px = si[idx - n] * precon[idx - n];
                    double py = li[idx - n] * precon[idx - n];
                    e -= (py * py + tau * px * py);
                }

                if (e < sigma * Adiag[idx]) {
                    e = Adiag[idx];
                }

                // don't encase a small amount of fluid cells within solids, it messes with the preconditioner. But only modifying it only if e is big enough patches up the problem
                if (e * e > 1e-18) {
                    precon[idx] = 1.0 / std::sqrt(e);
                }
            }
        }
    }

    void applyPreconditioner(std::vector<double> *dst, std::vector<double> *a) {
        for (int i = 1; i < numX - 1; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int32_t idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;

                double t = (*a)[idx];

                // if cell - S is fluid, then do: (si - S) * (precon - S) * (dst - S)
                if (fluid_attributes.cellType[idx - 1] == FLUID_CELL) {
                    t -= si[idx - 1] * precon[idx - 1] * (*dst)[idx - 1];
                }

                // if cell - L is fluid, then do: (si - L) * (precon - L) * (dst - L)
                if (fluid_attributes.cellType[idx - n] == FLUID_CELL) {
                    t -= li[idx - n] * precon[idx - n] * (*dst)[idx - n];
                }

                (*dst)[idx] = t * precon[idx];
            }
        }

        for (int i = numX - 2; i > 0; --i) {
            for (int j = numY - 2; j > 0; --j) {
                int32_t idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;

                double t = (*dst)[idx];

                // if cell + S is fluid, then do: (si) * (precon) * (dst + S)
                if (fluid_attributes.cellType[idx + 1] == FLUID_CELL) {
                    t -= si[idx] * precon[idx] * (*dst)[idx + 1];
                }

                // if cell + L is fluid, then do: (li) * (precon) * (dst + L)
                if (fluid_attributes.cellType[idx + n] == FLUID_CELL) {
                    t -= li[idx] * precon[idx] * (*dst)[idx + n];
                }

                (*dst)[idx] = t * precon[idx];
            }
        }
    }

    void matVec(std::vector<double> *dst, std::vector<double> *b) {
        for (int i = 1; i < numX - 1; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int32_t idx = i * n + j;
                
                double t = Adiag[idx] * (*b)[idx];

                // (si - S) * (b - S)
                t += si[idx - 1] * (*b)[idx - 1];

                // (li - L) * (b - L)
                t += li[idx - n] * (*b)[idx - n];

                // (si) * (b + S)
                t += si[idx] * (*b)[idx + 1];

                // (li) * (b + L)
                t += li[idx] * (*b)[idx + n];

                (*dst)[idx] = t;
            }
        }
    }

    void PCGproject() {
        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));

        std::fill(begin(pressure), end(pressure), 0.f);

        setUpResidual();
        setUpA();
        buildPreconditioner();
        applyPreconditioner(&z, &residual);
        std::copy(begin(z), end(z), begin(search));

        // DotMulti needs to be debugged, EqualsPlusTimesMulti good

        /*double sigma = 0.0;
        Dot(&z, &residual, 0, numX * numY, sigma);*/
        auto start = std::chrono::high_resolution_clock::now();

        double sigma = DotMulti(&z, &residual);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(DotTime, duration.count());

        for (int iter = 0; iter < numPressureIters && sigma > 0; ++iter) {
            matVec(&z, &search);

            //double denom = 0.0;
            //Dot(&z, &search, 0, numX * numY, denom);
            //double alpha = sigma / denom;

            double alpha = sigma / DotMulti(&z, &search);

            start = std::chrono::high_resolution_clock::now();
            ScaledAdd(pressure, search, alpha);
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            addValueToAverage(scaledAddTime, duration.count());

            ScaledAdd(residual, z, -alpha);

            start = std::chrono::high_resolution_clock::now();
            applyPreconditioner(&z, &residual);           // applying the preconditioner is literally the ONLY thing making this so slow
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            addValueToAverage(preconditioningTime, duration.count());

            //double sigmaNew = 0.0;
            //Dot(&z, &residual, 0, numX * numY, sigmaNew);

            double sigmaNew = DotMulti(&z, &residual);
            
            start = std::chrono::high_resolution_clock::now();
            EqualsPlusTimesMulti(&search, &z, sigmaNew / sigma);
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            addValueToAverage(EqualsPlusTime, duration.count());

            sigma = sigmaNew;
        }

        applyPressureMulti();
    }

    void SORproject() {
        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));

        std::fill(begin(pressure), end(pressure), 0.f);

        for (int iter = 0; iter < numPressureIters; ++iter) {
            for (int i = 1; i < this->numX - 1; ++i) {
                for (int j = 1; j < this->numY - 1; ++j) {
                    if (fluid_attributes.cellType[i * n + j] != FLUID_CELL) continue;

                    float leftType = fluid_attributes.cellType[(i - 1) * n + j] <= AIR_CELL ? 1 : 0;
                    float rightType = fluid_attributes.cellType[(i + 1) * n + j] <= AIR_CELL ? 1 : 0;
                    float topType = fluid_attributes.cellType[i * n + j - 1] <= AIR_CELL ? 1 : 0;
                    float bottomType = fluid_attributes.cellType[i * n + j + 1] <= AIR_CELL ? 1 : 0;

                    float divideBy = leftType + rightType + topType + bottomType;
                    if (divideBy == 0.f) continue;

                    float divergence = fluid_attributes.u[(i + 1) * n + j] - fluid_attributes.u[i * n + j] + fluid_attributes.v[i * n + j + 1] - fluid_attributes.v[i * n + j];

                    if (fluid_attributes.particleRestDensity > 0.f) {
                        float compression = fluid_attributes.cellDensities[i * n + j] - fluid_attributes.particleRestDensity;
                        if (compression > 0.f) {
                            divergence -= k * compression;
                        }
                    }

                    float p = divergence / divideBy;
                    p *= overRelaxation;

                    fluid_attributes.u[i * n + j] += leftType * p;
                    fluid_attributes.u[(i + 1) * n + j] -= rightType * p;
                    fluid_attributes.v[i * n + j] += topType * p;
                    fluid_attributes.v[i * n + j + 1] -= bottomType * p;
                }
            }
        }
    }

    float curl(int i, int j) {
        const int32_t n = this->numY;
        const float denom = 1.f / (2.f * fluid_attributes.cellSpacing);
        const int32_t leftType = fluid_attributes.cellType[(i - 1) * n + j] == FLUID_CELL;
        const int32_t rightType = fluid_attributes.cellType[(i + 1) * n + j] == FLUID_CELL;
        const int32_t topType = fluid_attributes.cellType[i * n + j - 1] == FLUID_CELL;
        const int32_t bottomType = fluid_attributes.cellType[i * n + j + 1] == FLUID_CELL;
        if (!leftType || !rightType || !topType || !bottomType) {
            return 0.f;
        }
        return ((fluid_attributes.v[(i + 1) * n + j] * bottomType - fluid_attributes.v[(i - 1) * n + j] * topType) - (fluid_attributes.u[i * n + j + 1] * rightType - fluid_attributes.u[i * n + j - 1] * leftType)) * denom;
    }

    float calcVorticity(int i, int j) {
        float curl = this->curl(i, j);
        return std::abs(curl);// * curl;
    }

    void calcVorticityConfinement(bool red, int32_t startColumn, int32_t endColumn) {
        const int32_t n = this->numY;
        for (int32_t i = startColumn; i < endColumn; ++i) {
            for (int32_t j = 1; j < numY - 1; ++j) {
                if (red) {
                    if ((i + j) % 2 != 0) continue; 
                }
                else {
                    if ((i + j) % 2 == 0) continue; 
                }
                float dx = abs(curl(i, j + 1)) - abs(curl(i, j - 1));
                float dy = abs(curl(i + 1, j)) - abs(curl(i - 1, j));

                const float len = std::sqrt(dx * dx + dy * dy);

                const float invLen = 1.f / (len + (len == 0.f)) - (len == 0.f);

                dx *= invLen;
                dy *= invLen;

                const float c = curl(i, j);

                fluid_attributes.v[i * n + j] += c * dx * dt * fluid_attributes.vorticityStrength;
                fluid_attributes.u[i * n + j] += c * dy * dt * fluid_attributes.vorticityStrength;
            }
        }
    }

    void applyVorticityConfinementRedBlack() {
        const int32_t numColumnsPerThread = (numX - 2) / numThreads;
        const int32_t numMissedColumns = numX - 2 - numColumnsPerThread * numThreads;

        for (int i = 0; i < numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->calcVorticityConfinement(true, 1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->calcVorticityConfinement(true, numX - 1 - numMissedColumns, numX - 1);

        fluid_attributes.thread_pool.waitForCompletion();

        for (int i = 0; i < numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->calcVorticityConfinement(false, 1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->calcVorticityConfinement(false, numX - 1 - numMissedColumns, numX - 1);

        fluid_attributes.thread_pool.waitForCompletion();
    }

    void includeRigidObject(const bool mouseDown, const bool justPressed) {
        const float extend = 20 * fluid_attributes.cellSpacing;
        if (mouseDown) {
            float vx = (objectX - objectPrevX) * 100;
            float vy = (objectY - objectPrevY) * 100;
            for (int i = 1; i < numX - 1; i++) {
                for (int j = 1; j < numY - 1; j++) {
                    int cellNr = i * n + j;
                    if (fluid_attributes.cellType[i * n + j] == SOLID_CELL) continue;
                    float dx = (i + 0.5) * fluid_attributes.cellSpacing - objectX;
                    float dy = (j + 0.5) * fluid_attributes.cellSpacing - objectY;

                    if (dx * dx + dy * dy < objectRadius * objectRadius + extend) {
                        //fluid_attributes.cellType[i * n + j] = FLUID_CELL;
                        
                        if (fluid_attributes.cellType[cellNr - n] != SOLID_CELL) {
                            fluid_attributes.u[cellNr] = vx;
                        }
                        if (fluid_attributes.cellType[cellNr + n] != SOLID_CELL) {
                            fluid_attributes.u[cellNr + n] = vx;
                        }
                        if (fluid_attributes.cellType[cellNr - 1] != SOLID_CELL) {
                            fluid_attributes.v[cellNr] = vy;
                        }
                        if (fluid_attributes.cellType[cellNr + 1] != SOLID_CELL) {
                            fluid_attributes.v[cellNr + 1] = vy;
                        }
                    }
                }
            }
            objectPrevX = objectX;
            objectPrevY = objectY;
            objectX = std::max(fluid_attributes.cellSpacing + objectRadius, std::min(mouseX, fluid_attributes.WIDTH - fluid_attributes.cellSpacing - objectRadius));
            objectY = std::max(fluid_attributes.cellSpacing + objectRadius, std::min(mouseY, fluid_attributes.HEIGHT - fluid_attributes.cellSpacing - objectRadius));
            if (justPressed) {
                objectPrevX = objectX;
                objectPrevY = objectY;
            }
        }
        else {
            objectX = std::max(fluid_attributes.cellSpacing + objectRadius, std::min(objectX, fluid_attributes.WIDTH - fluid_attributes.cellSpacing - objectRadius));
            objectY = std::max(fluid_attributes.cellSpacing + objectRadius, std::min(objectY, fluid_attributes.HEIGHT - fluid_attributes.cellSpacing - objectRadius));
            objectPrevX = objectX;
            objectPrevY = objectY;
        }
    }

    void drawRigidObject(sf::RenderWindow& window) {
        this->objectDrawer.setPosition(objectX, objectY);
        window.draw(this->objectDrawer);
    }

    void generate() {
        float separation = radius * 2.1;
        int wideNum = std::floor((2 * generatorRadius) / (separation));
        int highNum = wideNum;
        int numPotentiallyAdded = wideNum * highNum;

        float generatorRightBound = mouseX + generatorRadius;
        float generatorBottomBound = mouseY + generatorRadius;

        float starting_px = std::max(mouseX - generatorRadius + radius, fluid_attributes.cellSpacing + radius);
        float starting_py = std::max(mouseY - generatorRadius + radius, fluid_attributes.cellSpacing + radius);
        float px = starting_px;
        float py = starting_py;
        bool offset = true;

        fluid_attributes.positions.resize(2 * fluid_attributes.num_particles + 2 * numPotentiallyAdded);
        int addedTo = 0;

        for (int i = 0; i < numPotentiallyAdded; ++i) {
            float prevPx = px;
            float prevPy = py;
            if (prevPy > generatorBottomBound || prevPy + radius > fluid_attributes.HEIGHT - fluid_attributes.cellSpacing) {
                break;
            }

            int cellX = px / fluid_attributes.cellSpacing;
            int cellY = py / fluid_attributes.cellSpacing;
            int cellNr = cellX * n + cellY;

            px += separation;
            if (px > generatorRightBound) {
                px = starting_px;
                /*if (offset) {
                    px += 0.5 * separation;
                }*/
                py += separation;
                offset = !offset;
            }

            if (fluid_attributes.cellType[cellNr] != AIR_CELL || prevPx - radius < fluid_attributes.cellSpacing || prevPx + radius > fluid_attributes.WIDTH - fluid_attributes.cellSpacing || prevPy - radius < fluid_attributes.cellSpacing) {
                continue;
            }

            fluid_attributes.positions[2 * (fluid_attributes.num_particles + addedTo)] = prevPx;
            fluid_attributes.positions[2 * (fluid_attributes.num_particles + addedTo) + 1] = prevPy;

            addedTo++;
        }

        fluid_attributes.num_particles += addedTo;

        fluid_attributes.positions.resize(2 * fluid_attributes.num_particles);
        fluid_attributes.velocities.resize(2 * fluid_attributes.num_particles);
        fluid_renderer.particleColors.resize(3 * fluid_attributes.num_particles);
        fluid_renderer.va.resize(4 * fluid_attributes.num_particles);

        int start = fluid_attributes.num_particles - addedTo;

        for (int i = start; i < fluid_attributes.num_particles; i++) {
            int idx1 = 2 * i;
            int idx2 = 3 * i;
            int idx3 = 4 * i;
            int idx4 = i - start;

            fluid_attributes.velocities[idx1] = 0.f;
            fluid_attributes.velocities[idx1 + 1] = 0.f;

            fluid_renderer.particleColors[idx2] = 255;
            fluid_renderer.particleColors[idx2 + 1] = 255;
            fluid_renderer.particleColors[idx2 + 2] = 255;

            fluid_renderer.va[idx3].texCoords = {0.f, 0.f};
            fluid_renderer.va[idx3 + 1].texCoords = {fluid_renderer.texture_size.x, 0.f};
            fluid_renderer.va[idx3 + 2].texCoords = {fluid_renderer.texture_size.x, fluid_renderer.texture_size.y};
            fluid_renderer.va[idx3 + 3].texCoords = {0.f, fluid_renderer.texture_size.y};
        }

        fluid_attributes.particlesPerThread = fluid_attributes.num_particles / numThreads;
        fluid_attributes.numMissedParticles = fluid_attributes.num_particles - numThreads * particlesPerThread;

        this->collisions.resize(fluid_attributes.num_particles);
        fluid_attributes.temperatures.resize(fluid_attributes.num_particles);

        this->nr0.resize(2 * fluid_attributes.num_particles);
        this->nr1.resize(2 * fluid_attributes.num_particles);
        this->nr2.resize(2 * fluid_attributes.num_particles);
        this->nr3.resize(2 * fluid_attributes.num_particles);

        this->d0.resize(2 * fluid_attributes.num_particles);
        this->d1.resize(2 * fluid_attributes.num_particles);
        this->d2.resize(2 * fluid_attributes.num_particles);
        this->d3.resize(2 * fluid_attributes.num_particles);

        //std::cout << fluid_attributes.num_particles << "\n";
    }

    void remove() {
        const int32_t numCovered = std::ceil(generatorRadius / scalingFactor);
        const uint32_t mouseColumn = std::floor(mouseX / scalingFactor);
        const uint32_t mouseRow = std::floor(mouseY / scalingFactor);
        const size_t double_len = fluid_attributes.positions.size();
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
    
                    if (2 * particleIndex + 2 < double_len) {
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
    
        fluid_attributes.positions.resize(2 * fluid_attributes.num_particles);
        fluid_attributes.velocities.resize(2 * fluid_attributes.num_particles);
        fluid_renderer.particleColors.resize(3 * fluid_attributes.num_particles);
        fluid_renderer.vaCopy.resize(4 * fluid_attributes.num_particles);
    
        fluid_renderer.va.resize(4 * fluid_attributes.num_particles);
        for (int i = 0; i < vaSize; ++i) {
            fluid_renderer.va[i] = fluid_renderer.vaCopy[i];
        }
    
        fluid_attributes.particlesPerThread = fluid_attributes.num_particles / numThreads;
        fluid_attributes.numMissedParticles = fluid_attributes.num_particles - numThreads * particlesPerThread;
    
        this->collisions.resize(fluid_attributes.num_particles);
        fluid_attributes.temperatures.resize(fluid_attributes.num_particles);
    
        this->nr0.resize(2 * fluid_attributes.num_particles);
        this->nr1.resize(2 * fluid_attributes.num_particles);
        this->nr2.resize(2 * fluid_attributes.num_particles);
        this->nr3.resize(2 * fluid_attributes.num_particles);
    
        this->d0.resize(2 * fluid_attributes.num_particles);
        this->d1.resize(2 * fluid_attributes.num_particles);
        this->d2.resize(2 * fluid_attributes.num_particles);
        this->d3.resize(2 * fluid_attributes.num_particles);
    }


    void drawGenerator(sf::RenderWindow& window) {
        generatorDrawer.setPosition(mouseX, mouseY);
        window.draw(generatorDrawer);
    }

    void drawForceObject(sf::RenderWindow& window) {
        forceObjectDrawer.setPosition(mouseX, mouseY);
        window.draw(forceObjectDrawer); 
    }

    void drawPencil(sf::RenderWindow& window) {
        int localX = static_cast<int>(mouseX / fluid_attributes.cellSpacing);
        int localY = static_cast<int>(mouseY / fluid_attributes.cellSpacing);
        int idx = localX * n + localY;

        float drawPosX = localX * fluid_attributes.cellSpacing + halfSpacing;
        float drawPosY = localY * fluid_attributes.cellSpacing + halfSpacing;

        if (leftMouseDown || !rightMouseDown) {
            pencil.setFillColor(sf::Color(0, 150, 0));
        }
        else {
            pencil.setFillColor(sf::Color(150, 0, 0));
        }

        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                if (localX + i > 0 && localY + j > 0 && localX + i < numX - 1 && localY + j < numY - 1) {
                    pencil.setPosition(drawPosX + i * fluid_attributes.cellSpacing, drawPosY + j * fluid_attributes.cellSpacing);
                    window.draw(pencil);
                }
            }
        }
    }

    void DrawDivergences(sf::RenderWindow& window) {
        float maxDiv = 0.f;
        for (int i = 1; i < numX - 1; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int idx = i * numY + j;
                if (fluid_attributes.cellType[idx] == FLUID_CELL) {
                    float div = fabsf(fluid_attributes.u[(i + 1) * n + j] - fluid_attributes.u[idx] + fluid_attributes.v[idx + 1] - fluid_attributes.v[idx]);

                    if (div > maxDiv) {
                        maxDiv = div;
                    }
                }
            }
        }

        if (maxDiv == 0.f) maxDiv = 1e-5f;

        for (int i = 1; i < numX - 1; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;
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

    void drawObstacles(sf::RenderWindow& window) {
        window.draw(obstacleVa, obstacleStates);
    }

    void addValueToAverage(float& value, float newValue) {
        value += (newValue - value) / steps;
    }

    void update(float dt_, sf::RenderWindow& window, bool leftMouseDown, bool justPressed, bool rightMouseDown) {
        //auto start = std::chrono::high_resolution_clock::now();
        sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
        this->mouseX = mouse_pos.x;
        this->mouseY = mouse_pos.y;
        
        if (!fluid_attributes.stop || fluid_attributes.step) {
            this->simulate(dt_, leftMouseDown, justPressed, rightMouseDown);
            fluid_attributes.step = false;
        }

        this->render(window);

        /*auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        addValueToAverage(SimStepTime, duration.count());*/
    }

    void simulate(float dt_, bool leftMouseDown_, bool rightMouseDown_, bool justPressed) {
        ++steps;

        // order of need of implementation/optimization:
            // 1) refactor all this messy code
            // 2) incompressibility -- implement MGPCG
            // 2.5) make a sampleVelocity(point) function so that you can do RK2 advection easier
            // 3) make it so that you pass in static arrays instead of just numbers of particles in the main file
            // 4) level set & fast sweeping for separattion from obstacles
            // 5) implement implicit density projection
            // 6) move divergence view into a vertex array
            // 7) updateDensity -- Same idea as to grid -- low priority
        
        //auto start = std::chrono::high_resolution_clock::now();

        dt = dt_;

        this->integrateMulti();

        if (fluid_attributes.fireActive) {
            std::fill(begin(collisions), end(collisions), 0);
        }
        /*auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count());*/


        //start = std::chrono::high_resolution_clock::now();
        addObjectsToGrids();

        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(FillGridTime, duration.count());*/

        


        //start = std::chrono::high_resolution_clock::now();
        leftMouseDown = leftMouseDown_;
        rightMouseDown = rightMouseDown_;

        if (solidDrawing && leftMouseDown) {
            this->drawSolids();
        }
        else if (solidDrawing && rightMouseDown) {
            this->eraseSolids();
        }

        if (generatorActive && leftMouseDown) {
            this->generate();
        }
        else if (generatorActive && rightMouseDown) {
            this->remove();
        }
        
        if (forceObjectActive && leftMouseDown) {
            this->makeForceObjectQueries(-250); // pulling, -250
        }
        else if (forceObjectActive && rightMouseDown) { 
            this->makeForceObjectQueries(1000); // pushing, 1000
        }

        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count());*/


        


        //start = std::chrono::high_resolution_clock::now();
        solveCollisions();
        
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(CollisionTime, duration.count());*/





        //start = std::chrono::high_resolution_clock::now();

        this->collideSurfacesMulti();

        this->constrainWallsMulti();

        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ObstacleCollisionTime, duration.count());*/





        //start = std::chrono::high_resolution_clock::now();
        
        this->transfer(true);
        
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ToGridTime, duration.count());*/



        //start = std::chrono::high_resolution_clock::now();
        
        this->updateParticleDensity();
        
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(DensityUpdateTime, duration.count());*/



        //start = std::chrono::high_resolution_clock::now();
        
        if (rigidObjectActive) {
            this->includeRigidObject(leftMouseDown, justPressed);
        }


        if (fluid_attributes.vorticityStrength != 0) {
            this->applyVorticityConfinementRedBlack();
        }
        
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count());*/

        


        //start = std::chrono::high_resolution_clock::now();
        
        //this->PCGproject();
        this->SORproject();
        
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ProjectionTime, duration.count());*/




        //start = std::chrono::high_resolution_clock::now();
        
        this->transfer(false);
        
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ToParticlesTime, duration.count());*/
    }

    void render(sf::RenderWindow& window) {
        //auto start = std::chrono::high_resolution_clock::now();
        
        //this->drawCells(window);

        if (renderPattern == 0) {
            fluid_renderer.UpdateVaDiffusionMulti();
        }

        else if (renderPattern == 1) {
            fluid_renderer.UpdateVaVelocityMulti();
        }

        else if (renderPattern == 2) {
            fluid_renderer.UpdateVaVorticityMulti();
        }

        else if (renderPattern == 3) {
            fluid_renderer.UpdateVaTemperatureMulti();
        }
        else if (renderPattern == 4) {
            this->DrawDivergences(window);
        }

        this->drawObstacles(window);

        if (renderPattern != 4) {
            fluid_renderer.DrawParticles();
        }

        if (forceObjectActive) {
            this->drawForceObject(window);
        }
        else if (rigidObjectActive) {
            this->drawRigidObject(window);
        }
        else if (generatorActive) {
            this->drawGenerator(window);
        }

        if (solidDrawing) {
            this->drawPencil(window);
        }

        //this->drawUVGrids(window);
        /*auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(RenderingTime, duration.count());*/
    }

    float getCombinedTime() {
        return FillGridTime + miscellaneousTime + CollisionTime + ObstacleCollisionTime + ToGridTime + DensityUpdateTime + ProjectionTime + ToParticlesTime + RenderingTime;
    }

    float getSimStepTime() {
        return SimStepTime;
    }

    float getFillGridTime() {
        return FillGridTime;
    }

    float getMiscellaneousTime() {
        return miscellaneousTime;
    }

    float getCollisionTime() {
        return CollisionTime;
    }

    float getObstacleCollisionTime() {
        return ObstacleCollisionTime;
    }

    float getToGridTime() {
        return ToGridTime;
    }

    float getDensityUpdateTime() {
        return DensityUpdateTime;
    }

    float getProjectionTime() {
        return ProjectionTime;
    }

    float getToParticlesTime() {
        return ToParticlesTime;
    }

    float getRenderingTime() {
        return RenderingTime;
    }

    float getEqualsPlusTime() {
        return this->EqualsPlusTime;
    }

    float getDotTime() {
        return this->DotTime;
    }

    float getScaledAddTime() {
        return this->scaledAddTime;
    }

    float getPreconditionTime() {
        return this->preconditioningTime;
    }

    void addToForceObjectRadius(float add) {
        if (forceObjectRadius + add > 0) {
            this->forceObjectRadius += add;
            this->forceObjectDrawer.setOrigin(forceObjectRadius, forceObjectRadius);
            this->forceObjectDrawer.setRadius(forceObjectRadius);
            this->checkForceObjectseparationDist = (this->radius + forceObjectRadius) * (this->radius + forceObjectRadius);
        }
    }

    void addToGeneratorRadius(float add) {
        if (generatorRadius + add > 0) {
            this->generatorRadius += add;
            this->generatorDrawer.setOrigin(generatorRadius, generatorRadius);
            this->generatorDrawer.setSize(sf::Vector2f(2 * generatorRadius, 2 * generatorRadius));
        }
    }

    void addToGravityX(float add) {
        fluid_attributes.gravityX += add;
    }

    float getGravityX() {
        return fluid_attributes.gravityX;
    }

    void addToGravityY(float add) {
        fluid_attributes.gravityY += add;
    }

    float getGravityY() {
        return fluid_attributes.gravityY;
    }

    void addToDivergenceModifier(float add) {
        this->k += add;
    }

    float getDivergenceModifier() {
        return this->k;
    }

    void setNextRenderPattern() {
        this->renderPattern++;
        if (this->renderPattern > 4) {
            this->renderPattern = 0;
        }
    }

    void setForceObjectActive(bool active) {
        this->forceObjectActive = active;
    }

    void setGeneratorActive(bool active) {
        this->generatorActive = active;
    }

    void setRigidObjectActive(bool active) {
        this->rigidObjectActive = active;
    }

    void setFireActive(bool active) {
        this->fireActive = active;
    }

    float getFireActive() {
        return this->fireActive;
    }

    void addToNumPressureIters(int32_t add) {
        numPressureIters += add;
    }

    int32_t getNumPressureIters() {
        return numPressureIters;
    }

    int getNumX() {
        return this->numX;
    }

    int getNumY() {
        return this->numY;
    }

    void setSolidDrawer(bool set) {
        this->solidDrawing = set;
    }

    bool getPencilActive() {
        return this->solidDrawing;
    }

    int getPencilRadius() {
        return this->pencilRadius;
    }

    void addToPencilRadius(int add) {
        this->pencilRadius += add;
    }

    bool getForceObjectActive() {
        return this->forceObjectActive;
    }

    bool getGeneratorActive() {
        return this->generatorActive;
    }

    void drawUVGrids(sf::RenderWindow& window) {
        sf::VertexArray line(sf::Lines, 2);
        int32_t n = numY;
        for (int i = 0; i < numX; ++i) {
            for (int j = 0; j < numY; ++j) {

                // draw u lines (left right)
                float uX = fluid_attributes.cellSpacing + i * fluid_attributes.cellSpacing;
                float uY = 1.5 * fluid_attributes.cellSpacing + j * fluid_attributes.cellSpacing;
                line[0].position = sf::Vector2f(uX, uY);
                line[0].color  = sf::Color(255, 150, 0);
                line[1].position = sf::Vector2f(uX + fluid_attributes.u[i * n + j], uY);
                line[1].color = sf::Color(255, 150, 0);
                window.draw(line);

                //draw v lines (top bottom)
                float vX = 1.5 * fluid_attributes.cellSpacing + i * fluid_attributes.cellSpacing;
                float vY = fluid_attributes.cellSpacing + j * fluid_attributes.cellSpacing;
                line[0].position = sf::Vector2f(vX, vY);
                line[0].color  = sf::Color::Red;
                line[1].position = sf::Vector2f(vX, vY + fluid_attributes.v[i * n + j]);
                line[1].color = sf::Color::Red;
                window.draw(line);
            }
        }
    }
};
