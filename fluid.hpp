#pragma once
#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <array>

#include "thread_pool.hpp"
#include "collision_grid.hpp"

float dot(const std::vector<float>& a, const std::vector<float>& b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.f);
}

void axpy(float alpha, const std::vector<float>& x, std::vector<float>& y) {
    for (size_t i = 0; i < x.size(); ++i) {
        y[i] += alpha * x[i];
    }
}

void scal(float alpha, std::vector<float>& x) {
    for (float& val : x) {
        val *= alpha;
    }
}

class Fluid {
    int numX;
    int numY;
    float numCells;
    float invSpacing;
    float cellSpacing;
    int numParticles;
    float radius;
    std::vector<float> u;
    std::vector<float> v;
    std::vector<float> du;
    std::vector<float> dv;
    std::vector<float> prevU;
    std::vector<float> prevV;
    std::vector<float> p;
    std::vector<float> cellType;
    std::vector<float> cellColor;
    std::vector<float> positions;
    std::vector<float> velocities;
    std::vector<float> particleDensity;
    float particleRestDensity = 0;
    bool interacting = false;
    int WIDTH;
    int HEIGHT;
    std::vector<sf::Color> colors;

    static constexpr float restitution = 0.f;
    float checkSeperationDist;
    float moveDist;
    int U_FIELD = 0;
    int V_FIELD = 1;
   
    int FLUID_CELL = 0;
    int AIR_CELL = 1;
    int SOLID_CELL = 2;

    sf::CircleShape particleDrawer;

    float objectRadius;
    float objectX;
    float objectY;
    float objectPrevX;
    float objectPrevY;
    sf::CircleShape objectDrawer;
    std::map<int, float> rigidObjectQueries;
    float checkRigidObjectSeperationDist;
    float objectXVel = 0;
    float objectYVel = 0;

    float mouseX = 0;
    float mouseY = 0;

    float forceObjectRadius = 250; // 200
    std::map<int, float> forceObjectQueries;
    sf::CircleShape forceObjectDrawer;
    float checkForceObjectSeperationDist;

    float generatorRadius = forceObjectRadius;
    sf::RectangleShape generatorDrawer;

    std::vector<int> particleColors;

    float gravity;

    std::array<std::array<int, 3>, 100> gradient;
    std::array<std::array<int, 3>, 4> colorMap{{{0, 51, 102}, {0, 153, 204}, {102, 255, 204}, {255, 255, 255}}};

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
    // dark blue: {0, 0, 128}, {0, 128, 255}, {255, 128, 0}, {255, 255, 0}
    // lightning mcqueen: {255, 0, 0}, {255, 69, 0}, {255, 165, 0}, {255, 255, 0}
    // rainbow: {255, 0, 0}, {255, 255, 0}, {0, 255, 0}, {0, 200, 255} 

    sf::VertexArray va{sf::PrimitiveType::Quads};

    sf::Texture texture;

    sf::RenderStates states;

    float k;

    float timeForTransfer;
    float timeForSeperation;
    float timeForIncompressibility;

    int numRowsPerThread;
    int numMissedRows;
    
    tp::ThreadPool& thread_pool;

    const float colorDiffusionCoeff = 0.001f;

    float diffusionRatio;

    float scalingFactor;
    int32_t scaledWIDTH;
    int32_t scaledHEIGHT;

    CollisionGrid grid;

    float dt;

    bool forceObjectActive = true;
    bool rigidObjectActive = false;
    bool generatorActive = false;

    float textureSizeX;
    float textureSizeY;

    float seperationInit;

    float vorticityStrength;

    std::vector<float> confinementForce;

    std::vector<float> vorticity;
    std::vector<float> vorticityMag;

    // any variables or methods commented out are for constant memory spatial hasing 
    /*int tableSize;
    int numObjects;
    std::vector<int> cellCount;
    std::vector<int> particleArray;
    std::vector<int> allHashCells;
    float spacing;*/

    float flipRatio;
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

    std::vector<float> temperatures;
    float groundConductivity = 30.f; // 20
    float interConductivity = 15.f;  // 15
    float fireStrength = 100.f;      // 75
    float tempDiffusion = 0.1f;      // 0.1

    int32_t renderPattern = 0;

    bool fireActive = false;

public:
    Fluid(float WIDTH, float HEIGHT, float cellSpacing, int numParticles, float gravity, float k, float diffusionRatio, float seperationInit, float vorticityStrength_, float flipRatio_, float overRelaxation_, float numPressureIters_, tp::ThreadPool& tp)
        : numX(std::floor(WIDTH / cellSpacing)), numY(std::floor(HEIGHT / cellSpacing)), numCells(numX * numY), numParticles(numParticles), WIDTH(WIDTH), HEIGHT(HEIGHT), gravity(gravity), k(k), diffusionRatio(diffusionRatio), seperationInit(seperationInit), vorticityStrength(vorticityStrength_), flipRatio(flipRatio_), overRelaxation(overRelaxation_), numPressureIters(numPressureIters_), thread_pool(tp) {

            this->temperatures.resize(numParticles);

            this->nr0.resize(2 * numParticles);
            this->nr1.resize(2 * numParticles);
            this->nr2.resize(2 * numParticles);
            this->nr3.resize(2 * numParticles);

            this->d0.resize(2 * numParticles);
            this->d1.resize(2 * numParticles);
            this->d2.resize(2 * numParticles);
            this->d3.resize(2 * numParticles);

            this->u.resize(numCells);
            this->v.resize(numCells);
            this->du.resize(numCells);
            this->dv.resize(numCells);
            this->prevU.resize(numCells);
            this->prevV.resize(numCells);
            this->p.resize(numCells);
            this->cellType.resize(numCells);
            this->cellColor.resize(3 * numCells);
            this->positions.resize(2 * numParticles);
            this->velocities.resize(2 * numParticles);
            this->particleDensity.resize(numCells);
            this->particleColors.resize(3 * numParticles);
            std::fill(begin(particleColors), end(particleColors), 0);
            this->va.resize(numParticles * 4);
            texture.loadFromFile("white_circle.png");
            texture.generateMipmap();
            auto const texture_size = static_cast<sf::Vector2f>(texture.getSize());
            for (int index = 0; index < numParticles; ++index) {
                int i = 4 * index;
                va[i].texCoords = {0.f, 0.f};
                va[i + 1].texCoords = {texture_size.x, 0.f};
                va[i + 2].texCoords = {texture_size.x, texture_size.y};
                va[i + 3].texCoords = {0.f, texture_size.y};
            }
            states.texture = &texture;

            textureSizeX = texture_size.x;
            textureSizeY = texture_size.y;

            this->cellSpacing = std::max(WIDTH / numX, HEIGHT / numY);
            this->invSpacing = 1.f / this->cellSpacing;

            this->radius = 0.3 * cellSpacing;
            this->particleDrawer.setOutlineThickness(0.f);
            this->particleDrawer.setRadius(this->radius);
            this->particleDrawer.setOrigin(this->radius, this->radius);
            this->particleDrawer.setFillColor(sf::Color(0, 150, 255));

            this->scalingFactor = 2 * radius;

            this->scaledWIDTH = std::ceil(static_cast<float>(WIDTH) / scalingFactor);
            this->scaledHEIGHT = std::ceil(static_cast<float>(HEIGHT) / scalingFactor);

            grid = CollisionGrid(scaledWIDTH, scaledHEIGHT);

            // initializing particle positions
            float seperation = radius / seperationInit;  // gridsize: 65: 2.5; 70: 2.3; 90: 2.f; 130: 1.2

            int wideNum = std::floor((WIDTH - 2 * cellSpacing - 2) / (radius * seperation));
            int highNum = numParticles / wideNum;

            float starting_px = radius + cellSpacing + 2;//(WIDTH - (radius * seperation * wideNum)) / 2 + radius;
            float starting_py = (HEIGHT - (radius * seperation * highNum)) / 2 + radius;

            float px = starting_px;
            float py = starting_py;

            int addTo = numParticles - (wideNum * highNum);

            bool offset = true;
            for (int i = 0; i < wideNum * highNum + addTo; ++i) {
                this->positions[i * 2] = px;
                this->positions[i * 2 + 1] = py;
                this->particleColors[3 * i] = 0;
                this->particleColors[3 * i + 1] = 0;
                this->particleColors[3 * i + 2] = 255;

                px += this->radius * seperation;

                if ((i + 1) % wideNum == 0) {
                    px = starting_px;
                    if (offset) {
                        px += this->radius;
                    }
                    py += this->radius * seperation;
                    offset = !offset;
                }
            }

            int n = this->numY;
            for (int i = 0; i < this->numX; ++i) {
                for (int j = 0; j < this->numY; ++j) {
                    if (i == 0 || j == 0 || i == this->numX - 1 || j == this->numY - 1) {
                        this->cellType[i * n + j] = SOLID_CELL;
                    }
                }
            }

            this->moveDist = 2 * radius;
            this->checkSeperationDist = moveDist * moveDist;

            this->objectRadius = cellSpacing * 3;
            this->objectX = std::floor(WIDTH / 2);
            this->objectY = 2 * objectRadius + 10;
            this->objectPrevX = objectX;
            this->objectPrevY = objectY;
            this->objectDrawer.setOrigin(objectRadius, objectRadius);
            this->objectDrawer.setRadius(objectRadius);
            this->objectDrawer.setFillColor(sf::Color(255, 0, 0));
            this->checkRigidObjectSeperationDist = (this->radius + objectRadius) * (this->radius + objectRadius);

            this->forceObjectDrawer.setOrigin(forceObjectRadius, forceObjectRadius);
            this->forceObjectDrawer.setRadius(forceObjectRadius);
            this->forceObjectDrawer.setOutlineThickness(1.f);
            this->forceObjectDrawer.setFillColor(sf::Color::Transparent);
            this->forceObjectDrawer.setOutlineColor(sf::Color::Red); 
            this->checkForceObjectSeperationDist = (this->radius + forceObjectRadius) * (this->radius + forceObjectRadius);

            this->generatorDrawer.setOrigin(generatorRadius, generatorRadius);
            this->generatorDrawer.setSize(sf::Vector2f(2 * generatorRadius, 2 * generatorRadius));
            this->generatorDrawer.setOutlineThickness(1.f);
            this->generatorDrawer.setFillColor(sf::Color::Transparent);
            this->generatorDrawer.setOutlineColor(sf::Color::Red); 

            // linearly interpolate between the values in colorMap to create a gradient array 
            float num_colors = colorMap.size() - 1; // number of colors - 1
            float num_steps = 1.f * gradient.size() / num_colors; //num_steps = 50 * key_range
            int index = 0;
            for (int i = 0; i < num_colors; ++i) {  
                for (int x = 0; x < num_steps; ++x) {
                    float t = 1.f * x / num_steps;  // Interpolation factor
                    // Linear interpolation for r, g, b values between colorMap[i] andcolorMap [i+1]
                    int r = (int)(colorMap[i][0] * (1 - t) + colorMap[i + 1][0] * t);
                    int g = (int)(colorMap[i][1] * (1 - t) + colorMap[i + 1][1] * t);
                    int b = (int)(colorMap[i][2] * (1 - t) + colorMap[i + 1][2] * t);
                    gradient[index] = std::array<int, 3>{r, g, b};
                    index++;
                }
            }

            num_colors = tempMap.size() - 1; // number of colors - 1
            num_steps = 1.f * tempgradient.size() / num_colors; //num_steps = 50 * key_range
            index = 0;
            for (int i = 0; i < num_colors; ++i) {  
                for (int x = 0; x < num_steps; ++x) {
                    float t = 1.f * x / num_steps;  // Interpolation factor
                    // Linear interpolation for r, g, b values between colorMap[i] andcolorMap [i+1]
                    int r = (int)(tempMap[i][0] * (1 - t) + tempMap[i + 1][0] * t);
                    int g = (int)(tempMap[i][1] * (1 - t) + tempMap[i + 1][1] * t);
                    int b = (int)(tempMap[i][2] * (1 - t) + tempMap[i + 1][2] * t);
                    tempgradient[index] = std::array<int, 3>{r, g, b};
                    index++;
                }
            }

            /*this->spacing = 2 * this->radius;
            this->tableSize = 2 * numParticles;
            this->cellCount.resize(this->tableSize + 1);
            this->particleArray.resize(numParticles);
            this->allHashCells.resize(numParticles);*/
    }

    float clamp(float x, float min, float max) {
        if (x < min) {
            return min;
        }
        else if (x > max) {
            return max;
        }
        return x;
    }

    void createRandomPositions() {
        std::uniform_int_distribution<int> randWidth(radius, WIDTH - radius);
        std::uniform_int_distribution<int> randHeight(radius, HEIGHT - radius);

        std::random_device rd;
        std::mt19937 mt(rd());

        for (int i = 0; i < numParticles; ++i) {
            positions[i * 2] = randWidth(mt);
            positions[i * 2 + 1] = randHeight(mt);
        }
    }

    void solveIncompressibilityRedBlackForRows(const uint32_t startColumn, const uint32_t endColumn, const bool red) {
        int32_t n = this->numY;
        for (int i = startColumn; i < endColumn; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                if (red) { 
                    if ((i + j) % 2 != 0) {
                        continue;
                    }
                }
                else { 
                    if ((i + j) % 2 == 0) {
                        continue;
                    }
                } 
                float leftType = cellType[(i - 1) * n + j] <= AIR_CELL ? 1 : 0;
                float rightType = cellType[(i + 1) * n + j] <= AIR_CELL ? 1 : 0;
                float topType = cellType[i * n + j - 1] <= AIR_CELL ? 1 : 0;
                float bottomType = cellType[i * n + j + 1] <= AIR_CELL ? 1 : 0;
                float divideBy = leftType + rightType + topType + bottomType;

                if (divideBy == 0.f) continue;

                float divergence = this->u[(i + 1) * n + j] - this->u[i * n + j] + this->v[i * n + j + 1] - this->v[i * n + j];

                if (this->particleRestDensity > 0.f) {
                    float compression = this->particleDensity[i * n + j] - this->particleRestDensity;
                    if (compression > 0.f) {
                        divergence -= this->k * compression;
                    }
                }

                float p = divergence / divideBy;
                p *= overRelaxation;

                this->u[i * n + j] += leftType * p;
                this->u[(i + 1) * n + j] -= rightType * p;
                this->v[i * n + j] += topType * p;
                this->v[i * n + j + 1] -= bottomType * p;
            }
        }
    }

    /*void solveIncompressibilityMulti() {
        for (int32_t iter = 0; iter < )
        thread_pool.addTask([this, i, slice_size]{
            int32_t const start{2 * i * slice_size};
            int32_t const end  {start + slice_size};
            solveCollisionThreaded(start, end);
        });
    }*/

    void integrate(const uint32_t startIndex, const uint32_t endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            this->positions[2 * i] += this->velocities[2 * i] * dt;
            this->positions[2 * i + 1] += this->velocities[2 * i + 1] * dt;
            this->velocities[2 * i + 1] += gravity * dt;

            if (this->fireActive && this->renderPattern == 2 && positions[2 * i + 1] < HEIGHT - cellSpacing - 10) {
                this->velocities[2 * i + 1] -= fireStrength * temperatures[i] * dt;
            }
            if (this->temperatures[i] > 0) {
                this->temperatures[i] -= tempDiffusion;
            }
        }
    }

    void addObjectsToGrid()
    {
        grid.clear();

        const float minX = cellSpacing;
        const float maxX = WIDTH - cellSpacing;
        const float minY = cellSpacing;
        const float maxY = HEIGHT - cellSpacing;

        uint32_t i{0};
        for (int32_t index = 0; index < numParticles; ++index) {
            float x = positions[2 * index];
            float y = positions[2 * index + 1];
            if (x > minX && x < maxX && y > minY && y < maxY) {
                grid.addAtom(static_cast<int32_t>(x / scalingFactor), static_cast<int32_t>(y / scalingFactor), i);
            }
            ++i;
        }
    }

    void solveContact(uint32_t index, uint32_t otherIndex)
    {
        constexpr float eps           = 0.0001f;
        const float o2_o1X = positions[2 * index] - positions[2 * otherIndex];
        const float o2_o1Y = positions[2 * index + 1] - positions[2 * otherIndex + 1];

        const float dist2 = o2_o1X * o2_o1X + o2_o1Y * o2_o1Y;

        if (dist2 < checkSeperationDist && dist2 > eps) {
            const float dist          = sqrt(dist2);
            const float delta = 0.5f * (moveDist - dist) / dist;
            const float col_vecX = o2_o1X * delta;
            const float col_vecY = o2_o1Y * delta;

            positions[2 * index] += col_vecX;
            positions[2 * index + 1] += col_vecY;

            positions[2 * otherIndex] -= col_vecX;
            positions[2 * otherIndex + 1] -= col_vecY;

            const float transfer = (temperatures[index] - temperatures[otherIndex]) * interConductivity * dt;
            temperatures[index] -= transfer;
            temperatures[otherIndex] += transfer;
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
                checkAtomCellCollisions(atom_idx, grid.data[index + grid.height + side]);
            }
            for (int32_t side = 0; side < 2; ++side) {
                checkAtomCellCollisions(atom_idx, grid.data[index + side]);   
            }
            //checkAtomCellCollisions(atom_idx, grid.data[index - grid.height]);
            checkAtomCellCollisions(atom_idx, grid.data[index - grid.height + 1]);
        }
    }

    void solveCollisionThreaded(uint32_t start, uint32_t end)
    {
        for (uint32_t idx{start}; idx < end; ++idx) {
            processCell(grid.data[idx], idx);
        }
    }

    void makeForceObjectQueries(const int32_t strength) {
        const int32_t numCovered = std::ceil(forceObjectRadius / scalingFactor);

        const uint32_t mouseColumn = std::floor(mouseX / scalingFactor);
        const uint32_t mouseRow = std::floor(mouseY / scalingFactor);

        for (int32_t i = -numCovered; i < numCovered + 1; ++i) {
            for (int32_t j = -numCovered; j < numCovered + 1; ++j) {

                if (mouseRow + i <= 1 || mouseRow + i >= scaledHEIGHT - 1 || mouseColumn + j <= 1 || mouseColumn + j >= scaledWIDTH - 1) continue;

                const auto cell = grid.data[mouseRow + i + grid.height * (mouseColumn + j)];

                for (uint32_t i{0}; i < cell.objects_count; ++i) {
                    const uint32_t particleIndex = cell.objects[i];

                    float dx = positions[particleIndex * 2] - mouseX;
                    float dy = positions[particleIndex * 2 + 1] - mouseY;
                    float d2 = dx * dx + dy * dy;

                    if (d2 > checkForceObjectSeperationDist || d2 == 0.0f) continue;

                    float d = std::sqrt(d2);

                    float edgeT = d / forceObjectRadius;
            
                    float centerT = 1 - edgeT;

                    velocities[2 * particleIndex] += (dx * strength - velocities[2 * particleIndex]) * centerT * dt;
                    velocities[2 * particleIndex + 1] += (dy * strength - velocities[2 * particleIndex + 1]) * centerT * dt;
                }
            }
        }
    }

    // Find colliding atoms
    void solveCollisions()
    {
        // Multi-thread grid
        const uint32_t thread_count = thread_pool.m_thread_count;
        const uint32_t slice_count  = thread_count * 2;
        const uint32_t slice_size   = (grid.width / slice_count) * grid.height;
        const uint32_t last_cell    = 2 * thread_count * slice_size;
        // Find collisions in two passes to avoid data races

        // First collision pass
        for (uint32_t i{0}; i < thread_count; ++i) {
            thread_pool.addTask([this, i, slice_size]{
                uint32_t const start{2 * i * slice_size};
                uint32_t const end  {start + slice_size};
                solveCollisionThreaded(start, end);
            });
        }
        // Eventually process rest if the world is not divisible by the thread count
        if (last_cell < grid.data.size()) {
            thread_pool.addTask([this, last_cell]{
                solveCollisionThreaded(last_cell, static_cast<uint32_t>(grid.data.size()));
            });
        }
        thread_pool.waitForCompletion();
        // Second collision pass
        for (uint32_t i{0}; i < thread_count; ++i) {
            thread_pool.addTask([this, i, slice_size]{
                uint32_t const start{(2 * i + 1) * slice_size};
                uint32_t const end  {start + slice_size};
                solveCollisionThreaded(start, end);
            });
        }
        thread_pool.waitForCompletion();
    }

    void constrainWalls(const uint32_t startIndex, const uint32_t endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            if (this->positions[2 * i] - radius < this->cellSpacing) {
                this->positions[2 * i] = radius + this->cellSpacing;
                if (this->velocities[2 * i] < 0) {
                    this->velocities[2 * i] *= -restitution;
                }
            }
            else if (this->positions[2 * i] + radius > WIDTH - this->cellSpacing) {
                this->positions[2 * i] = WIDTH - radius - this->cellSpacing;
                if (this->velocities[2 * i] > 0) {
                    this->velocities[2 * i] *= -restitution;
                }
            }
            if (this->positions[2 * i + 1] - radius < this->cellSpacing) {
                this->positions[2 * i + 1] = radius + this->cellSpacing;
                if (this->velocities[2 * i + 1] < 0) {
                    this->velocities[2 * i + 1] *= -restitution;
                }
            }
            else if (this->positions[2 * i + 1] + radius > HEIGHT - this->cellSpacing) {
                this->positions[2 * i + 1] = HEIGHT - radius - this->cellSpacing;
                if (this->velocities[2 * i + 1] > 0) {
                    this->velocities[2 * i + 1] *= -restitution;
                }
                const float remove = 5;
                if (this->temperatures[i] < tempgradient.size() && this->positions[2 * i] > WIDTH / remove && this->positions[2 * i] < WIDTH - WIDTH / remove) {
                    this->temperatures[i] += groundConductivity;
                }
            }
        }
    }

    void cacheTransferNodes2(int32_t start, int32_t end, float invHeight, int32_t component) {
        const float n = this->numY;
        const float h2 = invHeight;

        float dx = (component != 0) * h2;
        float dy = (component == 0) * h2;

        for (int32_t i = start; i < end; ++i) {
            float x = this->positions[2 * i];
            float y = this->positions[2 * i + 1];
            x = this->clamp(x, cellSpacing, (this->numX - 1) * cellSpacing);
            y = this->clamp(y, cellSpacing, (this->numY - 1) * cellSpacing);
            // x0 is the grid position to the left of the particle, x1 is the position to the right of  the particle. Both can only go up to the second to last cell to the right in the     gridbecause we dont want to be changing wall velocities
            int x0 = std::max(1, std::min((int)(std::floor((x - dx) * invSpacing)), this->numX - 1));
            // basically x - xCell to get the weight of that cell in relation to the particle
            // in this case, x is moved over to x - dx, and xCell is just grid position of x multiplied     by grid spacing
            float tx = ((x - dx) - x0 * cellSpacing) * invSpacing;
            // add 1 to get the cell to the right
            int x1 = std::min(x0 + 1, this->numX - 2);
            // this fixes a bug that makes water touching the left wall and ceiling explode sometimes 
            if (component == 0 && x0 == 1) {
                x0 = x1;
            }
            if (component == 1 && x0 == 1) {
                x1 = x0;
            }
            // same thing with y
            int y0 = std::max(0, std::min((int)(std::floor((y - dy) * invSpacing)), this->numY - 2));
            float ty = ((y - dy) - y0 * cellSpacing) * invSpacing;
            int y1 = std::min(y0 + 1, this->numY - 1);
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
        float h2 = 0.5 * cellSpacing;

        const int32_t numThreads = thread_pool.m_thread_count;
        const int32_t numParticlesPerThread = numParticles / numThreads;
        const int32_t numMissedParticles = numParticles - numThreads * numParticlesPerThread;

        for (int32_t i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i](){
                cacheTransferNodes2(numParticlesPerThread * i, numParticlesPerThread * i + numParticlesPerThread, h2, 0);
                cacheTransferNodes2(numParticlesPerThread * i, numParticlesPerThread * i + numParticlesPerThread, h2, 1);
            });
        }

        thread_pool.waitForCompletion();
    }

    void transferMulti(const bool& toGrid) {
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
        const float h2 = 0.5 * cellSpacing;

            // make a copy of the grid
        std::copy(std::begin(this->u), std::end(this->u), std::begin(this->prevU));
        std::copy(std::begin(this->v), std::end(this->v), std::begin(this->prevV));
        std::fill(begin(this->du), end(this->du), 0.f);
        std::fill(begin(this->dv), end(this->dv), 0.f);
        std::fill(begin(this->u), end(this->u), 0.f);
        std::fill(begin(this->v), end(this->v), 0.f);

            // initialize outside cells to solid, every inside cell to air
        for (int i = 0; i < this->numX; ++i) {
            for (int j = 0; j < this->numY; ++j) {
                    /*if (this->cellType[i * n + j] == SOLID_CELL) {
                        this->cellColor[3 * (i * n + j)] = 100;
                        this->cellColor[3 * (i * n + j) + 1] = 100;
                        this->cellColor[3 * (i * n + j) + 2] = 100;
                    }*/
                if (this->cellType[i * n + j] != SOLID_CELL) {
                        this->cellType[i * n + j] = AIR_CELL;
                       
                        /*this->cellColor[3 * (i * n + j)] = 0;
                        this->cellColor[3 * (i * n + j) + 1] = 250;
                        this->cellColor[3 * (i * n + j) + 2] = 255;
                        */
                }
            }
        }

            // initialize all cells that particles are in to fluid
        for (int i = 0; i < this->numParticles; ++i) {
            float x = this->positions[2 * i];
            float y = this->positions[2 * i + 1];

                // get cell coords
            int xi = this->clamp(std::floor(x * invSpacing), 0, this->numX - 1);
            int yi = this->clamp(std::floor(y * invSpacing), 0, this->numY - 1);

            int cellNr = xi * n + yi;
                // if a cell has particle(s) in it, change it to fluid cell
            if (this->cellType[cellNr] == AIR_CELL) {
                this->cellType[cellNr] = FLUID_CELL;
                   
                    /*
                    this->cellColor[3 * cellNr] = 0;
                    this->cellColor[3 * cellNr + 1] = 150;
                    this->cellColor[3 * cellNr + 2] = 255;
                    */
            }
        }
    }

    void transferToParticle(int32_t startIndex, int32_t endIndex) {
        const int32_t n = this->numY;
        for (int32_t i = startIndex; i < endIndex; ++i) {
            const int32_t nr0_u = nr0[2 * i];
            const int32_t nr1_u = nr1[2 * i];
            const int32_t nr2_u = nr2[2 * i];
            const int32_t nr3_u = nr3[2 * i];

            const float d0_u = d0[2 * i];
            const float d1_u = d1[2 * i];
            const float d2_u = d2[2 * i];
            const float d3_u = d3[2 * i];

            const int32_t nr0_v = nr0[2 * i + 1];
            const int32_t nr1_v = nr1[2 * i + 1];
            const int32_t nr2_v = nr2[2 * i + 1];
            const int32_t nr3_v = nr3[2 * i + 1];

            const float d0_v = d0[2 * i + 1];
            const float d1_v = d1[2 * i + 1];
            const float d2_v = d2[2 * i + 1];
            const float d3_v = d3[2 * i + 1];

            const float pvx = this->velocities[2 * i];
            const float pvy = this->velocities[2 * i + 1];
           
                // these will be used to make sure that air cells are not considered when transferring velocities back to particles
                // nr0 - offset is the same thing as [(i-1) * n]
                // if u is being considered, then we only have to check left and right cells ([nr0] and [nr0 - n])
                // if v is being considered, then we only have to check above and below cells ([nr0] and [nr0 - 1])
            const float valid0u = this->cellType[nr0_u] != AIR_CELL || this->cellType[nr0_u - n] != AIR_CELL;
            const float valid1u = this->cellType[nr1_u] != AIR_CELL || this->cellType[nr1_u - n] != AIR_CELL;
            const float valid2u = this->cellType[nr2_u] != AIR_CELL || this->cellType[nr2_u - n] != AIR_CELL;
            const float valid3u = this->cellType[nr3_u] != AIR_CELL || this->cellType[nr3_u - n] != AIR_CELL;

            const float valid0v = this->cellType[nr0_v] != AIR_CELL || this->cellType[nr0_v - 1] != AIR_CELL;
            const float valid1v = this->cellType[nr1_v] != AIR_CELL || this->cellType[nr1_v - 1] != AIR_CELL;
            const float valid2v = this->cellType[nr2_v] != AIR_CELL || this->cellType[nr2_v - 1] != AIR_CELL;
            const float valid3v = this->cellType[nr3_v] != AIR_CELL || this->cellType[nr3_v - 1] != AIR_CELL;

            const float divX = valid0u * d0_u + valid1u * d1_u + valid2u * d2_u + valid3u * d3_u;
            const float divY = valid0v * d0_v + valid1v * d1_v + valid2v * d2_v + valid3v * d3_v;

            float picV;
            float corr;
            float flipV;

            if (divX > 0.f) {
                picV = (valid0u * d0_u * this->u[nr0_u] + valid1u * d1_u * this->u[nr1_u] + valid2u * d2_u * this->u[nr2_u] + valid3u * d3_u * this->u[nr3_u]) / divX;

                corr = (valid0u * d0_u * (this->u[nr0_u] - this->prevU[nr0_u]) + valid1u * d1_u * (this->u[nr1_u] - this->prevU[nr1_u]) + valid2u * d2_u * (this->u[nr2_u] - this->prevU[nr2_u]) + valid3u * d3_u * (this->u[nr3_u] - this->prevU[nr3_u])) / divX;
                flipV = pvx + corr;
                this->velocities[2 * i] = (1.f - flipRatio) * picV + flipRatio * flipV;
            }

            if (divY > 0.f) {
                picV = (valid0v * d0_v * this->v[nr0_v] + valid1v * d1_v * this->v[nr1_v] + valid2v * d2_v * this->v[nr2_v] + valid3v * d3_v * this->v[nr3_v]) / divY;

                corr = (valid0v * d0_v * (this->v[nr0_v] - this->prevV[nr0_v]) + valid1v * d1_v * (this->v[nr1_v] - this->prevV[nr1_v]) + valid2v * d2_v * (this->v[nr2_v] - this->prevV[nr2_v]) + valid3v * d3_v * (this->v[nr3_v] - this->prevV[nr3_v])) / divY;
                flipV = pvy + corr;
                this->velocities[2 * i + 1] = (1.f - flipRatio) * picV + flipRatio * flipV;
            }
        }
    }

    void transferToParticles() {

        const int32_t numThreads = thread_pool.m_thread_count;
        const int32_t numParticlesPerThread = numParticles / numThreads;
        const int32_t numMissedParticles = numParticles - numParticlesPerThread * numThreads;

        // u and v grids are staggered, so make sure that you subtract half cell spacing from particle y positions when transferring this->u to particles and vice versa for this->v

        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->transferToParticle(i * numParticlesPerThread, i * numParticlesPerThread + numParticlesPerThread);
            });
        }

        this->transferToParticle(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();

        //this->transferToParticle(0, numParticles);
    }

    void transferToUGrid() {
        for (int32_t i = 0; i < this->numParticles; ++i) {
            const int32_t nr0_u = nr0[2 * i];
            const int32_t nr1_u = nr1[2 * i];
            const int32_t nr2_u = nr2[2 * i];
            const int32_t nr3_u = nr3[2 * i];

            const float d0_u = d0[2 * i];
            const float d1_u = d1[2 * i];
            const float d2_u = d2[2 * i];
            const float d3_u = d3[2 * i];

            const float pvx = this->velocities[2 * i];
            const float pvy = this->velocities[2 * i + 1];
           
            this->u[nr0_u] += pvx * d0_u;  
            this->u[nr1_u] += pvx * d1_u;
            this->u[nr2_u] += pvx * d2_u;
            this->u[nr3_u] += pvx * d3_u;

            this->du[nr0_u] += d0_u;
            this->du[nr1_u] += d1_u;
            this->du[nr2_u] += d2_u;
            this->du[nr3_u] += d3_u;
        }
    }

    void transferToVGrid() {
        for (int32_t i = 0; i < numParticles; ++i) {
            const int32_t nr0_v = nr0[2 * i + 1];
            const int32_t nr1_v = nr1[2 * i + 1];
            const int32_t nr2_v = nr2[2 * i + 1];
            const int32_t nr3_v = nr3[2 * i + 1];

            const float d0_v = d0[2 * i + 1];
            const float d1_v = d1[2 * i + 1];
            const float d2_v = d2[2 * i + 1];
            const float d3_v = d3[2 * i + 1];

            const float pvx = this->velocities[2 * i];
            const float pvy = this->velocities[2 * i + 1];
    
            this->v[nr0_v] += pvy * d0_v;  
            this->v[nr1_v] += pvy * d1_v;
            this->v[nr2_v] += pvy * d2_v;
            this->v[nr3_v] += pvy * d3_v;

            this->dv[nr0_v] += d0_v;
            this->dv[nr1_v] += d1_v;
            this->dv[nr2_v] += d2_v;
            this->dv[nr3_v] += d3_v;
        }
    }

    void transferToGrid() {
        const float n = this->numY;
        const float h2 = 0.5 * cellSpacing;

        // u and v grids are staggered, so make sure that you subtract half cell spacing from particle y positions when transferring this->u to particles and vice versa for this->v

        thread_pool.addTask([&]() {
            this->transferToUGrid();
        });
        thread_pool.addTask([&]() {
            this->transferToVGrid();
        });
        thread_pool.waitForCompletion();


        for (int i = 0; i < this->u.size(); ++i) {
            float prevNode = this->du[i];
            if (prevNode > 0.f) {
                this->u[i] /= prevNode;
            }
        }
        for (int i = 0; i < this->v.size(); ++i) {
            float prevNode = this->dv[i];
            if (prevNode > 0.f) {
                this->v[i] /= prevNode;
            }
        }
    }

    void updateParticleDensity() {
        int n = this->numY;
        float h = this->cellSpacing;
        float h1 = this->invSpacing;
        float h2 = 0.5 * h;

        std::fill(begin(this->particleDensity), end(this->particleDensity), 0);

        for (int i = 0; i < numParticles; ++i) {
            float x = this->positions[2 * i];
            float y = this->positions[2 * i + 1];

            x = this->clamp(x, h, (this->numX - 1) * h);
            y = this->clamp(y, h, (this->numY - 1) * h);
            // x0 is the grid position to the left of the particle, x1 is the position to the right of the particle. Both can only go up to the second to last cell to the right in the gridbecause we dont want to be changing wall velocities
            int x0 = std::max(1, std::min((int)(std::floor((x - h2) * h1)), this->numX - 2));
            // basically x - xCell to get the weight of that cell in relation to the particle
            // in this case, x is moved over to x - dx, and xCell is just grid position of x multiplied by grid spacing
            float tx = ((x - h2) - x0 * h) * h1;
            // add 1 to get the cell to the right
            int x1 = std::min(x0 + 1, this->numX - 1);
            int y0 = std::max(1, std::min((int)(std::floor((y - h2) * h1)), this->numY - 2));
            float ty = ((y - h2) - y0 * h) * h1;
            int y1 = std::min(y0 + 1, this->numY - 1);
           
            float sx = 1.f - tx;
            float sy = 1.f - ty;

            if (x0 < this->numX && y0 < this->numY) {
                this->particleDensity[x0 * n + y0] += sx * sy;
            }
            if (x1 < this->numX && y0 < this->numY) {
                this->particleDensity[x1 * n + y0] += tx * sy;
            }
            if (x1 < this->numX && y1 < this->numY) {
                this->particleDensity[x1 * n + y1] += tx * ty;
            }
            if (x0 < this->numX && y1 < this->numY) {
                this->particleDensity[x0 * n + y1] += sx * ty;
            }
        }

        if (this->particleRestDensity == 0.f) {
            float sum = 0.f;
            int numFluidCells = 0;

            for (int i = 0; i < this->numCells; ++i) {
                if (this->cellType[i] == FLUID_CELL) {
                    sum += this->particleDensity[i];
                    numFluidCells++;
                }
            }

            if (numFluidCells > 0) {
                this->particleRestDensity = sum / numFluidCells;
            }
        }
    }

    void solveIncompressibility() {
        std::fill(begin(this->p), end(this->p), 0);

        std::copy(std::begin(this->u), std::end(this->u), std::begin(this->prevU));
        std::copy(std::begin(this->v), std::end(this->v), std::begin(this->prevV));

        const int32_t n = this->numY;

        for (int32_t iter = 0; iter < numPressureIters; ++iter) {
            for (int32_t i = 1; i < this->numX - 1; ++i) {
                for (int32_t j = 1; j < this->numY - 1; ++j) {
                    if (this->cellType[i * n + j] != FLUID_CELL) continue;
                    // <= AIR_CELL just means "is either fluid or air cell" look at the constants
                    float leftType = cellType[(i - 1) * n + j] <= AIR_CELL ? 1 : 0;
                    float rightType = cellType[(i + 1) * n + j] <= AIR_CELL ? 1 : 0;
                    float topType = cellType[i * n + j - 1] <= AIR_CELL ? 1 : 0;
                    float bottomType = cellType[i * n + j + 1] <= AIR_CELL ? 1 : 0;

                    float divideBy = leftType + rightType + topType + bottomType;

                    if (divideBy == 0.f) continue;

                    float divergence = this->u[(i + 1) * n + j] - this->u[i * n + j] + this->v[i * n + j + 1] - this->v[i * n + j];

                    if (this->particleRestDensity > 0.f) {
                        float compression = this->particleDensity[i * n + j] - this->particleRestDensity;
                        if (compression > 0.f) {
                            divergence -= this->k * compression;
                        }
                    }

                    float p = divergence / divideBy;
                    p *= overRelaxation;

                    this->u[i * n + j] += leftType * p;
                    this->u[(i + 1) * n + j] -= rightType * p;
                    this->v[i * n + j] += topType * p;
                    this->v[i * n + j + 1] -= bottomType * p;
                }
            }
        }
    }

    float curl(int i, int j) {
        const int32_t n = this->numY;
        const int32_t leftType = cellType[(i - 2) * n + j] == FLUID_CELL;
        const int32_t rightType = cellType[(i + 1) * n + j] == FLUID_CELL;
        const int32_t topType = cellType[i * n + j - 2] == FLUID_CELL;
        const int32_t bottomType = cellType[i * n + j + 1] == FLUID_CELL;
        return (this->v[i * n + j + 1] * bottomType - this->v[i * n + j - 1] * topType) * 0.5f - (this->u[(i + 1) * n + j] * rightType - this->u[(i - 1) * n + j] * leftType) * 0.5f;
    }

    void calcVorticityConfinement(bool red, int32_t startColumn, int32_t endColumn) {
        const int32_t n = this->numY;
        for (int32_t i = startColumn; i < endColumn; ++i) {
            for (int32_t j = 1; j < numY - 1; ++j) {
                if (red) {
                    if ((i + j) % 2 != 0) continue; // Skip black cells
                }
                else {
                    if ((i + j) % 2 == 0) continue; // Skip red cells
                }
                float dx = abs(curl(i, j - 1)) - abs(curl(i, j + 1));
                float dy = abs(curl(i + 1, j)) - abs(curl(i - 1, j));

                const float len = std::sqrt(dx * dx + dy * dy);

                if (len > 0.f) {
                    dx /= len;
                    dy /= len; 
                }
                else {
                    dx = 0.f;
                    dy = 0.f;
                }

                const float c = curl(i, j);

                this->v[i * n + j] += c * dx * dt * vorticityStrength;
                this->u[i * n + j] += c * dy * dt * vorticityStrength;
            }
        }
    }

    void applyVorticityConfinementRedBlack() {
        const int32_t numThreads = thread_pool.m_thread_count;
        const int32_t numColumnsPerThread = (numX - 2) / numThreads;
        const int32_t numMissedColumns = numX - 2 - numColumnsPerThread * numThreads;

        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->calcVorticityConfinement(true, 1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->calcVorticityConfinement(true, numX - 1 - numMissedColumns, numX - 1);

        thread_pool.waitForCompletion();

        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->calcVorticityConfinement(false, 1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->calcVorticityConfinement(false, numX - 1 - numMissedColumns, numX - 1);

        thread_pool.waitForCompletion();
    }

    void includeRigidObject(const bool mouseDown, const bool justPressed) {
        if (mouseDown) {
            int n = numY;
            float vx = (objectX - objectPrevX) * 100;
            float vy = (objectY - objectPrevY) * 100;
            for (int i = 1; i < numX - 1; i++) {
                for (int j = 1; j < numY - 1; j++) {
                    //cellType[i * n + j] = AIR_CELL;
                    float dx = (i + 0.5) * cellSpacing - objectX;
                    float dy = (j + 0.5) * cellSpacing - objectY;

                    if (dx * dx + dy * dy < objectRadius *  objectRadius) {
                        cellType[i * n + j] = FLUID_CELL;
                        this->u[i*n + j] = vx;
                        this->u[(i+1)*n + j] = vx;
                        this->v[i*n + j] = vy;
                        this->v[i*n + j+1] = vy;
                    }
                }
            }
            objectPrevX = objectX;
            objectPrevY = objectY;
            objectX = std::max(cellSpacing + objectRadius, std::min(mouseX, WIDTH - cellSpacing - objectRadius));
            objectY = std::max(cellSpacing + objectRadius, std::min(mouseY, HEIGHT - cellSpacing - objectRadius));
            if (justPressed) {
                objectPrevX = objectX;
                objectPrevY = objectY;
            }
        }
        else {
            objectX = std::max(cellSpacing + objectRadius, std::min(objectX, WIDTH - cellSpacing - objectRadius));
            objectY = std::max(cellSpacing + objectRadius, std::min(objectY, HEIGHT - cellSpacing - objectRadius));
            objectPrevX = objectX;
            objectPrevY = objectY;
        }
    }

    void drawRigidObject(sf::RenderWindow& window) {
        this->objectDrawer.setPosition(objectX, objectY);
        window.draw(this->objectDrawer);
    }

    void generate(int32_t numParts) {
        this->positions.resize(numParticles + numParts);
        this->velocities.resize(numParticles + numParts);
        this->va.resize(4 * (numParticles + numParts));
        this->particleColors.resize(3 * (numParticles + numParts));
        for (int index = numParticles; index < numParticles + numParts; ++index) {
            int i = 4 * index;
            va[i].texCoords = {0.f, 0.f};
            va[i + 1].texCoords = {textureSizeX, 0.f};
            va[i + 2].texCoords = {textureSizeX, textureSizeY};
            va[i + 3].texCoords = {0.f, textureSizeY};
        }

        int wideNum = std::floor((2 * generatorRadius) / (2 * radius));
        int highNum = numParts / wideNum;
        float seperation = radius / 2.36f;
        int starting_px = mouseX - generatorRadius + radius;
        int starting_py = mouseY - generatorRadius + radius;
        int px = starting_px;
        int py = starting_py;
        int addTo = numParts - (wideNum * highNum);
        bool offset = true;
        for (int i = 0; i < wideNum * highNum + addTo; ++i) {
            this->positions[(i + numParticles) * 2] = px;
            this->positions[(i + numParticles) * 2 + 1] = py;
            this->velocities[(i + numParticles) * 2] = 0.f;
            this->velocities[(i + numParticles) * 2 + 1] = 0.f;
            px += this->radius * seperation;
            if ((i + 1) % wideNum == 0) {
                px = starting_px;
                if (offset) {
                    px += this->radius;
                }
                py += this->radius * seperation;
                offset = !offset;
            }
        }

        numParticles += numParts;
    }

    void drawGenerator(sf::RenderWindow& window) {
        generatorDrawer.setPosition(mouseX, mouseY);
        window.draw(generatorDrawer);
    }

    void applyForceObjectForces(float strength, float dt) {
        for (auto [otherParticleID, dist] : forceObjectQueries) {
            float dx = mouseX - positions[2 * otherParticleID];
            float dy = mouseY - positions[2 * otherParticleID + 1];
            float edgeT = dist / forceObjectRadius;
            float centerT = 1 - edgeT;

            velocities[2 * otherParticleID] += (dx * strength - velocities[2 * otherParticleID]) * centerT * dt;
            velocities[2 * otherParticleID + 1] += (dy * strength - velocities[2 * otherParticleID + 1]) * centerT * dt;
        } 
    }

    void drawForceObject(sf::RenderWindow& window) {
        forceObjectDrawer.setPosition(mouseX, mouseY);
        window.draw(forceObjectDrawer); 
        forceObjectQueries.clear();
    }

    void updateVertexArrayVelocity(uint32_t startIndex, uint32_t endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            int i = 4 * index;
            const float px = positions[2 * index];
            const float py = positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            sf::Color color;

            int vel = (int)(velocities[2 * index] * velocities[2 *  index] + velocities[2 * index + 1] * velocities[2 * index    + 1]) / 15000; 
            
            color = sf::Color(gradient[std::min(gradient.size() - 1, static_cast<unsigned long long>(vel))][0], gradient[std::min(gradient.size() - 1, static_cast<unsigned long long>(vel))][1], gradient[std::min(gradient.size() - 1, static_cast<unsigned long long>(vel))][2], 255);

            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void updateVertexArrayDiffusion(const uint32_t startIndex, const uint32_t endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            const float s = 0.01f;

            int i = 4 * index;
            const float px = positions[2 * index];
            const float py = positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            particleColors[3 * index] = clamp(this->particleColors[3    * index] - s, 0, 255);
            particleColors[3 * index + 1] = clamp(this->particleColors  [3 * index + 1] - s, 0, 255);
            particleColors[3 * index + 2] = clamp(this->particleColors  [3 * index + 2] + s, 0, 255);

            const int xi = clamp(std::floor(px * invSpacing), 1,    this->numX - 1);
            const int yi = clamp(std::floor(py * invSpacing), 1,    this->numY - 1);
            const int cellNr = xi * this->numY + yi;

            const float d0 = this->particleRestDensity;

            if (d0 > 0.f) {
                const float relDensity = this->particleDensity  [cellNr] / d0;
                if (relDensity < diffusionRatio) { // 1.25 for 20k w/   75 grid size, 1 for 30k w/75 grid size
                    particleColors[3 * index] = 204;
                    particleColors[3 * index + 1] = 204;
                    particleColors[3 * index + 2] = 255;
                }
            }

            sf::Color color = sf::Color(particleColors[3 * index],  particleColors[3 * index + 1], particleColors[3 * index +    2]);

            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void updateVertexArrayTemperature(const uint32_t startIndex, const uint32_t endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            const float s = 0.01f;

            int i = 4 * index;
            const float px = positions[2 * index];
            const float py = positions[2 * index + 1];

            float temp = temperatures[index];

            if (temp < 30 and fireActive) {
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

    void drawParticlesVertex(sf::RenderWindow& window) {
        window.draw(va, states);
    }

    void simulate(float dt_, sf::RenderWindow& window, bool leftMouseDown, bool justPressed, bool rightMouseDown) {

        // to grid, incompressibility
        dt = dt_;

        sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
        this->mouseX = mouse_pos.x;
        this->mouseY = mouse_pos.y;

        const bool mouseDown = leftMouseDown || rightMouseDown;

        if (generatorActive && mouseDown) {
            if (leftMouseDown) {
                this->generate(1);
            }
            /*else {
                this->remove();
            }*/
        }

        const uint32_t numThreads = thread_pool.m_thread_count;
        const uint32_t particlesPerThread = numParticles / numThreads;
        const uint32_t numMissedParticles = numParticles - numThreads * particlesPerThread;

        /*for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->integrate(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
            });
        }

        this->integrate(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();*/

        this->integrate(0, numParticles);

        //this->initializeSHConstantMem();
        addObjectsToGrid();

        // don't worry about this messing up the positions of the particles in relation to the grid; this only modifies their velocities
        if (forceObjectActive && mouseDown) {
            if (leftMouseDown) {
                this->makeForceObjectQueries(-250); // pulling, -250
            }
            else { 
                this->makeForceObjectQueries(1000); // pushing, 1000
            }
        }

        //auto start = std::chrono::high_resolution_clock::now();
        solveCollisions();
        /*auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "collision: " << duration.count() << " milliseconds" << "\n";*/

        constrainWalls(0, numParticles);

        /*for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->constrainWalls(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
            });
        }

        this->constrainWalls(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();*/

        //this->makeParticleQueriesConstantMem(0, numParticles);

        //this->makeForceObjectQueriesConstantMem(forceObjectActive);

        //auto start = std::chrono::high_resolution_clock::now();
        this->transferMulti(true);
        /*auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "to grid: " << duration.count() << " milliseconds" << "\n";*/

        this->updateParticleDensity();
        if (rigidObjectActive) {
            this->includeRigidObject(leftMouseDown, justPressed);
        }
        
        //this->solveIncompressibilityRedBlack(numPressureIters, overRelaxation);
        //this->applyVorticityConfinementRedBlack();
        
        //start = std::chrono::high_resolution_clock::now();
        this->solveIncompressibility();
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "incompressibility: " << duration.count() << " milliseconds" << "\n";*/

        // should be calculating vorticity confinement before the grid is solved, but calculating it after gives an artistic effect that I like, even if it's less physically accurate (vortex confinement isn't physically accurate anyways)
        //this->applyVorticityConfinementRedBlack();

        //start = std::chrono::high_resolution_clock::now();
        this->transferMulti(false);
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "to particles: " << duration.count() << " milliseconds" << "\n";*/

        //start = std::chrono::high_resolution_clock::now();
        if (renderPattern == 0) {
            for (int i = 0; i < numThreads; ++i) {
                thread_pool.addTask([&, i]() {
                    this->updateVertexArrayDiffusion(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
                });
            }

            this->updateVertexArrayDiffusion(numParticles - numMissedParticles, numParticles);

            thread_pool.waitForCompletion();
        }

        else if (renderPattern == 1) {
            for (int i = 0; i < numThreads; ++i) {
                thread_pool.addTask([&, i]() {
                    this->updateVertexArrayVelocity(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
                });
            }

            this->updateVertexArrayVelocity(numParticles - numMissedParticles, numParticles);

            thread_pool.waitForCompletion();
        }

        else if (renderPattern == 2) {
            for (int i = 0; i < numThreads; ++i) {
                thread_pool.addTask([&, i]() {
                    this->updateVertexArrayTemperature(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
                });
            }

            this->updateVertexArrayTemperature(numParticles - numMissedParticles, numParticles);

            thread_pool.waitForCompletion();
        }

        this->drawParticlesVertex(window);

        if (forceObjectActive) {
            this->drawForceObject(window);
        }
        else if (rigidObjectActive) {
            this->drawRigidObject(window);
        }
        else if (generatorActive) {
            this->drawGenerator(window);
        }
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "rendering: " << duration.count() << " milliseconds" << "\n";*/
    }

    void addToForceObjectRadius(float add) {
        if (forceObjectRadius + add > 0) {
            this->forceObjectRadius += add;
            this->forceObjectDrawer.setOrigin(forceObjectRadius, forceObjectRadius);
            this->forceObjectDrawer.setRadius(forceObjectRadius);
            this->checkForceObjectSeperationDist = (this->radius + forceObjectRadius) * (this->radius + forceObjectRadius);
        }
    }

    void addToGeneratorRadius(float add) {
        if (generatorRadius + add > 0) {
            this->generatorRadius += add;
            this->generatorDrawer.setOrigin(generatorRadius, generatorRadius);
            this->generatorDrawer.setSize(sf::Vector2f(2 * generatorRadius, 2 * generatorRadius));
        }
    }

    void addToGravity(float add) {
        this->gravity += add;
    }

    float getGravity() {
        return this->gravity;
    }

    void addToDivergenceModifier(float add) {
        this->k += add;
    }

    float getDivergenceModifier() {
        return this->k;
    }

    void addToVorticityStrength(float add) {
        this->vorticityStrength += add;
    }

    float getVorticityStrength() {
        return this->vorticityStrength;
    }

    float getTimeForSeperation() {
        return this->timeForSeperation;
    }

    float getTimeForInc() {
        return this->timeForIncompressibility;
    }

    float getTimeForTrans() {
        return this->timeForTransfer;
    }

    void setNextRenderPattern() {
        this->renderPattern++;
        if (this->renderPattern > 2) {
            this->renderPattern = 0;
        }
    }

    float getFlipRatio() {
        return this->flipRatio;
    }

    void addToFlipRatio(const float add) {
        this->flipRatio += add;
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

    /*void initializeSHConstantMem() {

        std::fill(this->cellCount.begin(), this->cellCount.end(), 0);
        std::fill(this->allHashCells.begin(), this->allHashCells.end(), 0);
        std::fill(this->particleArray.begin(), this->particleArray.end(), 0);
        // initialize cells in cellCount
        for (int i = 0; i < this->numParticles; ++i) {
            int hashedCell = hashCoords(this->positions[i * 2], this->positions[i * 2 + 1]);
            this->allHashCells[i] = hashedCell;
            this->cellCount[hashedCell]++;
        }
       
        // calc partial sum
        int sum = 0;
        for (int i = 0; i < this->tableSize; ++i) {
            sum += this->cellCount[i];
            this->cellCount[i] = sum;
        }
        this->cellCount[this->tableSize] = sum;
       
        // fill particle array
        for (int i = 0; i < this->numParticles; ++i) {
            int hashedCell = this->allHashCells[i];
            cellCount[hashedCell]--;
            particleArray[cellCount[hashedCell]] = i;
        }
    }

    void makeParticleQueriesConstantMem(int startIndex, int endIndex) {
        for (int c = startIndex; c < endIndex; ++c) {
            for (int i = -1; i < 2; ++i) {
                for (int j = -1; j < 2; ++j) {
                    int hashedCell = this->hashCoords(this->positions[c * 2] + i * this->spacing, this->positions[c * 2 + 1] + j * this->spacing * 0.5);
                    int start = this->cellCount[hashedCell];
                    int end = this->cellCount[hashedCell + 1];
                    for (int p = start; p < end; ++p) {
                        int otherParticleID = this->particleArray[p];
                        if (otherParticleID == c) continue;

                        float dx = this->positions[otherParticleID * 2] - this->positions[c * 2];
                        float dy = this->positions[otherParticleID * 2 + 1] - this->positions[c * 2 + 1];
                        float d2 = dx * dx + dy * dy;
                        if (d2 > checkSeperationDist || d2 == 0.0) continue;
                        float d = std::sqrt(d2); 
                        float s = 0.5 * (moveDist - d) / d;
                        dx *= s;
                        dy *= s;
                        this->positions[2 * c] -= dx;
                        this->positions[2 * c + 1] -= dy;
                        this->positions[2 * otherParticleID] += dx;
                        this->positions[2 * otherParticleID + 1] += dy;
                    }
                }
            }
        }
    }

    void makeForceObjectQueriesConstantMem(bool forceObjectActive) {
        if (forceObjectActive) {
            // might have to make this std::max(1, ...); 
            // using cellSpacing instead of spacing just works for the constant memory hashing
            float spacing = cellSpacing * 0.5;
            float numCovered = int(std::ceil(forceObjectRadius / spacing));
            for (int i = -numCovered; i < numCovered + 1; ++i) {
                for (int j = -numCovered; j < numCovered + 1; ++j) {
                    int hashedCell = this->hashCoords(mouseX + i * spacing, mouseY + j * spacing);
                    int start = this->cellCount[hashedCell];
                    int end = this->cellCount[hashedCell + 1];
                    for (int p = start; p < end; ++p) {
                        int otherParticleID = this->particleArray[p];
                        float dx = this->positions[otherParticleID * 2] - mouseX;
                        float dy = this->positions[otherParticleID * 2 + 1] - mouseY;
                        float d2 = dx * dx + dy * dy;
                        if (d2 > checkForceObjectSeperationDist || d2 == 0.0) continue;
                        float d = std::sqrt(d2);
                        forceObjectQueries[otherParticleID] = d;
                    }
                }
            }
        }
    }

    // would need a gigantic grid for incompressibility multithreading to be viable 
    void solveIncompressibilityRedBlack(const int numIters, const float overRelaxation) {
        std::fill(begin(this->p), end(this->p), 0);

        std::copy(std::begin(this->u), std::end(this->u), std::begin(this->prevU));
        std::copy(std::begin(this->v), std::end(this->v), std::begin(this->prevV));


        for (int iter = 0; iter < numIters; ++iter) {
            this->solvePass(overRelaxation, true);

            this->solvePass(overRelaxation, false);
        }
    }
    
    void drawHashNums(sf::RenderWindow& window, sf::Text& text) {
        for (int i = 0; i < WIDTH; i += this->spacing) {
            for (int j = 0; j < HEIGHT; j += this->spacing) {
                text.setPosition(i, j);
                text.setString(std::to_string(this->hashCoords(i, j)));
                window.draw(text);
            }
        }
    }
    
    
    /*int hashCoords(int xi, int yi) {
        int hash = (intCoord(xi) * 92837111) ^ (intCoord(yi) * 689287499);
        return std::abs(hash) % this->tableSize;
    }

    int intCoord(int coord) {
        // add 1 so that you dont get 0 for intCoord(xi) or intCoord(yi) in the hashCoords
        return std::floor(coord / this->spacing) + 1;
    }*/

    // Red-black gauss seidel implementation of solving incompressibility. 
    // May or may not converge faster than normal iteration, I just wanted to see if it would be viable for multithreading. 
    /*void solvePass(const float overRelaxation, const bool red) {
        // First Pass: Red Cells (i + j) % 2 == 0
        int n = numY;
        for (int i = 1; i < this->numX - 1; ++i) {
            for (int j = 1; j < this->numY - 1; ++j) {
                if (red) {
                    if ((i + j) % 2 != 0) continue; // Skip black cells
                }
                else {
                    if ((i + j) % 2 == 0) continue; // Skip red cells
                }
                if (this->cellType[i * n + j] != FLUID_CELL) continue;

                float leftType = cellType[(i - 1) * n + j] <= AIR_CELL ? 1 : 0;
                float rightType = cellType[(i + 1) * n + j] <= AIR_CELL ? 1 : 0;
                float topType = cellType[i * n + j - 1] <= AIR_CELL ? 1 : 0;
                float bottomType = cellType[i * n + j + 1] <= AIR_CELL ? 1 : 0;

                float divideBy = leftType + rightType + topType + bottomType;
                if (divideBy == 0.f) continue;

                float divergence = this->u[(i + 1) * n + j] - this->u[i * n + j] + this->v[i * n + j + 1] - this->v[i * n + j];

                if (this->particleRestDensity > 0.f) {
                    float k = 6.f;
                    float compression = this->particleDensity[i * n + j] - this->particleRestDensity;
                    if (compression > 0.f) {
                        divergence -= k * compression;
                    }
                }

                float p = divergence / divideBy;
                p *= overRelaxation;

                this->u[i * n + j] += leftType * p;
                this->u[(i + 1) * n + j] -= rightType * p;
                this->v[i * n + j] += topType * p;
                this->v[i * n + j + 1] -= bottomType * p;
            }
        }
    }*/
   
};
