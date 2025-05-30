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

    bool forceObjectActive = true;
    bool rigidObjectActive = false;
    bool generatorActive = false;

    float textureSizeX;
    float textureSizeY;

    float obstacleTextureSizeX;
    float obstacleTextureSizeY;

    int32_t renderPattern = 0;

    std::vector<uint32_t> collisions;

    sf::RectangleShape cellDrawer;

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

    PressureSolver &pressure_solver;

    FluidState &fluid_attributes;

    FluidRenderer &fluid_renderer;

    TransferGrid &transfer_grid;

    sf::CircleShape circleDrawer;

public:
    FluidHandler(float k, float overRelaxation_, float numPressureIters_, FluidState& fas, PressureSolver& ps, TransferGrid& tg, FluidRenderer& fr): k(k), fluid_attributes(fas), pressure_solver(ps), transfer_grid(tg), fluid_renderer(fr) {

            /*font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");
            text.setFont(font);
            text.setPosition(10, 10);
            text.setFillColor(sf::Color::White);*/

            circleDrawer.setRadius(fluid_attributes.cellSpacing / 6);    

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

            auto size = sf::Vector2f(fluid_attributes.cellSpacing, fluid_attributes.cellSpacing);
            this->cellDrawer.setSize(size);
            this->cellDrawer.setOutlineColor(sf::Color::White);
            this->cellDrawer.setOutlineThickness(1.f);

            this->collisions.resize(fluid_attributes.num_particles);

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
        }
    }

    void integrateMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->integrate(start, end);
        });
    }

    void makeFire(int start, int end) {
        for (int i = start; i < end; ++i) {
            if (fluid_attributes.positions[2 * i + 1] < fluid_attributes.HEIGHT - fluid_attributes.cellSpacing - 10) {
                fluid_attributes.velocities[2 * i + 1] -= fluid_attributes.fireStrength * fluid_attributes.temperatures[i] * dt;
                if (fluid_attributes.temperatures[i] > 0) {
                    fluid_attributes.temperatures[i] -= fluid_attributes.tempDiffusion * dt;
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
        const float maxX = fluid_attributes.WIDTH - fluid_attributes.cellSpacing;
        const float minY = fluid_attributes.cellSpacing;
        const float maxY = fluid_attributes.HEIGHT - fluid_attributes.cellSpacing;

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
            fluid_attributes.temperatures[index] -= transfer * dt;
            fluid_attributes.temperatures[otherIndex] += transfer * dt;

            collisions[index]++;
            collisions[otherIndex]++;


        }
    }

    void checkAtomCellCollisions(uint32_t atom_idx, const CollisionCell& c)
    {
        for (uint32_t i = 0; i < c.objects_count; ++i) {
            solveContact(atom_idx, c.objects[i]);
        }
    }

    void processCell(const CollisionCell& c, uint32_t index)
    {
        for (uint32_t i = 0; i < c.objects_count; ++i) {
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
        for (uint32_t idx = start; idx < end; ++idx) {
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
        const uint32_t slice_count  = fluid_attributes.numThreads * 2;
        const uint32_t slice_size   = (collisionGrid.width / slice_count) * collisionGrid.height;
        const uint32_t last_cell    = 2 * fluid_attributes.numThreads * slice_size;
        
        for (uint32_t i = 0; i < fluid_attributes.numThreads; ++i) {
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
        
        for (uint32_t i = 0; i < fluid_attributes.numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([this, i, slice_size]{
                uint32_t const start{(2 * i + 1) * slice_size};
                uint32_t const end  {start + slice_size};
                solveCollisionThreaded(start, end);
            });
        }
        fluid_attributes.thread_pool.waitForCompletion();
    }

    void heatGround(int32_t start, int32_t end) {
        for (int i = start; i < end; ++i) {
            if (fluid_attributes.positions[2 * i + 1] + radius > fluid_attributes.HEIGHT - fluid_attributes.cellSpacing) {
                const float remove = 10.f; // % of the floor from the sides that you want not heated
                if (renderPattern == 3 && fluid_attributes.temperatures[i] < fluid_renderer.tempgradient.size() && fluid_attributes.positions[2 * i] > fluid_attributes.WIDTH / remove && fluid_attributes.positions[2 * i] < fluid_attributes.WIDTH - fluid_attributes.WIDTH / remove) {
                    fluid_attributes.temperatures[i] += fluid_attributes.groundConductivity * dt;
                }
            }
        }
    }

    void heatGroundMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->heatGround(start, end);
        });
    }

    void constrainWalls(const uint32_t startIndex, const uint32_t endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            int xi = 2 * i;
            int yi = xi + 1;
            if (fluid_attributes.positions[xi] - radius < fluid_attributes.cellSpacing) {
                fluid_attributes.positions[xi] = radius + fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[xi] < 0) {
                    fluid_attributes.velocities[xi] = 0.f;
                }
            }
            else if (fluid_attributes.positions[xi] + radius > fluid_attributes.WIDTH - fluid_attributes.cellSpacing) {
                fluid_attributes.positions[xi] = fluid_attributes.WIDTH - radius - fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[xi] > 0) {
                    fluid_attributes.velocities[xi] = 0.f;
                }
            }
            if (fluid_attributes.positions[yi] - radius < fluid_attributes.cellSpacing) {
                fluid_attributes.positions[yi] = radius + fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[yi] < 0) {
                    fluid_attributes.velocities[yi] = 0.f;
                }
            }
            else if (fluid_attributes.positions[yi] + radius > fluid_attributes.HEIGHT - fluid_attributes.cellSpacing) {
                fluid_attributes.positions[yi] = fluid_attributes.HEIGHT - radius - fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[yi] > 0) {
                    fluid_attributes.velocities[yi] = 0.f;
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
        const int32_t numColumnsPerThread = (numX - 2) / fluid_attributes.numThreads;
        const int32_t numMissedColumns = numX - 2 - numColumnsPerThread * fluid_attributes.numThreads;

        for (int i = 0; i < fluid_attributes.numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->calcVorticityConfinement(true, 1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->calcVorticityConfinement(true, numX - 1 - numMissedColumns, numX - 1);

        fluid_attributes.thread_pool.waitForCompletion();

        for (int i = 0; i < fluid_attributes.numThreads; ++i) {
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

        fluid_attributes.positions.resize(2 * (fluid_attributes.num_particles + numPotentiallyAdded));
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

        fluid_attributes.particlesPerThread = fluid_attributes.num_particles / fluid_attributes.numThreads;
        fluid_attributes.numMissedParticles = fluid_attributes.num_particles - fluid_attributes.numThreads * fluid_attributes.particlesPerThread;

        this->collisions.resize(fluid_attributes.num_particles);
        fluid_attributes.temperatures.resize(fluid_attributes.num_particles);

        transfer_grid.nr0.resize(2 * fluid_attributes.num_particles);
        transfer_grid.nr1.resize(2 * fluid_attributes.num_particles);
        transfer_grid.nr2.resize(2 * fluid_attributes.num_particles);
        transfer_grid.nr3.resize(2 * fluid_attributes.num_particles);

        transfer_grid.d0.resize(2 * fluid_attributes.num_particles);
        transfer_grid.d1.resize(2 * fluid_attributes.num_particles);
        transfer_grid.d2.resize(2 * fluid_attributes.num_particles);
        transfer_grid.d3.resize(2 * fluid_attributes.num_particles);
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
        for (int i = 0; i < 4 * fluid_attributes.num_particles; ++i) {
            fluid_renderer.va[i] = fluid_renderer.vaCopy[i];
        }
    
        fluid_attributes.particlesPerThread = fluid_attributes.num_particles / fluid_attributes.numThreads;
        fluid_attributes.numMissedParticles = fluid_attributes.num_particles - fluid_attributes.numThreads * fluid_attributes.particlesPerThread;
    
        this->collisions.resize(fluid_attributes.num_particles);
        fluid_attributes.temperatures.resize(fluid_attributes.num_particles);
    
        transfer_grid.nr0.resize(2 * fluid_attributes.num_particles);
        transfer_grid.nr1.resize(2 * fluid_attributes.num_particles);
        transfer_grid.nr2.resize(2 * fluid_attributes.num_particles);
        transfer_grid.nr3.resize(2 * fluid_attributes.num_particles);
    
        transfer_grid.d0.resize(2 * fluid_attributes.num_particles);
        transfer_grid.d1.resize(2 * fluid_attributes.num_particles);
        transfer_grid.d2.resize(2 * fluid_attributes.num_particles);
        transfer_grid.d3.resize(2 * fluid_attributes.num_particles);
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

    void drawActiveUVNodes(sf::RenderWindow& window) {
        float dx = fluid_attributes.cellSpacing;

        for (int i = 0; i < numX; ++i) {
            for (int j = 0; j < numY; ++j) {
                int idx = i * n + j;

                bool activeSolidU = false;
                bool activeSolidV = false;
                bool airU = false;
                bool airV = false;

                if (i > 0 && ((fluid_attributes.cellType[idx] == FLUID_CELL && fluid_attributes.cellType[idx - n] == SOLID_CELL) || (fluid_attributes.cellType[idx - n] == FLUID_CELL && fluid_attributes.cellType[idx] == SOLID_CELL)))
                    activeSolidU = true;
                        
                if ((fluid_attributes.cellType[idx] == FLUID_CELL && fluid_attributes.cellType[idx - 1] == SOLID_CELL) || (fluid_attributes.cellType[idx - 1] == FLUID_CELL && fluid_attributes.cellType[idx] == SOLID_CELL))
                    activeSolidV = true;

                if (fluid_attributes.cellType[idx] != FLUID_CELL && fluid_attributes.cellType[idx - n] != FLUID_CELL) 
                    airU = true;

                if (fluid_attributes.cellType[idx] != FLUID_CELL && fluid_attributes.cellType[idx - 1] != FLUID_CELL)
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

    void drawObstacles(sf::RenderWindow& window) {
        window.draw(obstacleVa, obstacleStates);
    }

    void addValueToAverage(float& value, float newValue) {
        value += (newValue - value) / steps;
    }

    void update(float dt_, sf::RenderWindow& window, bool leftMouseDown, bool justPressed, bool rightMouseDown) {
        sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
        this->mouseX = mouse_pos.x;
        this->mouseY = mouse_pos.y;
        
        if (!fluid_attributes.stop || fluid_attributes.step) {
            //auto start = std::chrono::high_resolution_clock::now();
            
            this->simulate(dt_, leftMouseDown, justPressed, rightMouseDown);
            fluid_attributes.step = false;

            /*auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            addValueToAverage(SimStepTime, duration.count());*/
        }

        this->render(window);

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

        if (renderPattern == 3 && fluid_attributes.fireActive) {
            this->makeFireMulti();
        }

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

        if (renderPattern == 3) {
            this->heatGroundMulti();
        }

        this->collideSurfacesMulti();

        this->constrainWallsMulti();

        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ObstacleCollisionTime, duration.count());*/





        auto start = std::chrono::high_resolution_clock::now();
    
        transfer_grid.TransferToGrid();
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ToGridTime, duration.count());



        //start = std::chrono::high_resolution_clock::now();
        
        transfer_grid.updateParticleDensity();
        
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
        
        //pressure_solver.projectMICCG(1);
        //pressure_solver.projectSOR(pressure_solver.numPressureIters);
        //pressure_solver.projectRedBlackGS(pressure_solver.numPressureIters);
        pressure_solver.projectRedBlackGSMulti(pressure_solver.numPressureIters, fluid_attributes.numThreads);
        //pressure_solver.projectCG(pressure_solver.numPressureIters);

        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ProjectionTime, duration.count());*/




        //start = std::chrono::high_resolution_clock::now();
        
        transfer_grid.TransferToParticles();
        
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ToParticlesTime, duration.count());*/
    }

    void render(sf::RenderWindow& window) {
        auto start = std::chrono::high_resolution_clock::now();
        
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
            this->drawActiveUVNodes(window);
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

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(RenderingTime, duration.count());
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

    void setNextRenderPattern() {
        this->renderPattern++;
        if (this->renderPattern > 4) {
            this->renderPattern = 0;
        }
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

    void setRigidObjectActive(bool active) {
        this->rigidObjectActive = active;
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

};
