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

    static constexpr float restitution = 0.f;
    float checkseparationDist;
    float moveDist;
   
    int FLUID_CELL = 0;
    int AIR_CELL = 1;
    int SOLID_CELL = 2;

    float objectRadius;
    float objectX;
    float objectY;
    float objectPrevX;
    float objectPrevY;
    sf::CircleShape objectDrawer;
    std::map<int, float> rigidObjectQueries;
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

    std::vector<int> particleColors;

    float gravityX;
    float gravityY;

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

    sf::VertexArray va{sf::PrimitiveType::Quads};

    sf::Texture texture;

    sf::RenderStates states;

    float k;

    float timeForTransfer;
    float timeForseparation;
    float timeForIncompressibility;

    int numRowsPerThread;
    int numMissedRows;
    
    tp::ThreadPool& thread_pool;

    const float colorDiffusionCoeff = 0.001f;

    float diffusionRatio;

    float scalingFactor;
    int32_t scaledWIDTH;
    int32_t scaledHEIGHT;

    CollisionGrid collisionGrid;
    CollisionGrid cellOccupantsGrid;

    float dt;

    bool forceObjectActive = true;
    bool rigidObjectActive = false;
    bool generatorActive = false;

    float textureSizeX;
    float textureSizeY;

    // see if radius / 2.36 works
    float separationInit;

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
    // how quickly the ground transfers heat to the particles
    float groundConductivity = 100.f; // 30
    // how quickly particles transfer heat between themselves
    float interConductivity = 25.f;  // 15
    // how quickly particles accelerate upwards due to heat
    float fireStrength = 75.f;      // 100
    // how quickly heated particles lose that heat
    float tempDiffusion = 0.1f;      // 0.1

    int32_t renderPattern = 0;

    bool fireActive = false;

    std::vector<uint32_t> collisions;

    std::vector<double> direction;
    std::vector<double> residual;
    std::vector<double> pressure;
    std::vector<double> Ad;

    std::vector<double> Adiag;
    std::vector<double> si;
    std::vector<double> li;
    std::vector<double> precon;
    std::vector<double> z;
    std::vector<double> search;

    sf::RectangleShape cellDrawer;

    bool stop = false;
    bool step = false;

    uint32_t numThreads;
    uint32_t particlesPerThread;
    uint32_t numMissedParticles; 

    bool solidDrawing = false;

    int n;

    sf::Font font;

    sf::Text text;

    bool leftMouseDown = false;
    bool rightMouseDown = false;

    float halfSpacing;

    //sf::CircleShape circleDrawer;

    std::vector<sf::Vertex> vaCopy;

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

public:
    Fluid(float WIDTH, float HEIGHT, float cellSpacing, int numParticles, float gravityX_, float gravityY_, float k, float diffusionRatio, float separationInit, float vorticityStrength_, float flipRatio_, float overRelaxation_, float numPressureIters_, tp::ThreadPool& tp)
        : numX(std::floor(WIDTH / cellSpacing)), numY(std::floor(HEIGHT / cellSpacing)), numCells(numX * numY), numParticles(numParticles), WIDTH(WIDTH), HEIGHT(HEIGHT), gravityX(gravityX_), gravityY(gravityY_), k(k), diffusionRatio(diffusionRatio), separationInit(separationInit), vorticityStrength(vorticityStrength_), flipRatio(flipRatio_), overRelaxation(overRelaxation_), numPressureIters(numPressureIters_), thread_pool(tp) {

            /*circleDrawer.setFillColor(sf::Color(255, 0, 0));
            circleDrawer.setRadius(5);
            circleDrawer.setOrigin(5.f / 2.f, 5.f / 2.f);*/

            /*font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");
            text.setFont(font);
            text.setPosition(10, 10);
            text.setFillColor(sf::Color::White);*/

            this->cellSpacing = std::max(WIDTH / numX, HEIGHT / numY);
            this->invSpacing = 1.f / this->cellSpacing;

            this->radius = 0.3 * cellSpacing;

            this->halfSpacing = cellSpacing / 2;

            pencil.setSize(sf::Vector2f{cellSpacing, cellSpacing});
            pencil.setOrigin(halfSpacing, halfSpacing);
            pencil.setOutlineThickness(1);
            pencil.setOutlineColor(sf::Color::Black);

            this->numThreads = thread_pool.m_thread_count;
            this->particlesPerThread = numParticles / numThreads;
            this->numMissedParticles = numParticles - numThreads * particlesPerThread;

            this->dotProducts.resize(numThreads);

            auto size = sf::Vector2f(cellSpacing, cellSpacing);
            this->cellDrawer.setSize(size);
            this->cellDrawer.setOutlineColor(sf::Color::White);
            this->cellDrawer.setOutlineThickness(1.f);

            this->Adiag.resize(numX * numY);
            this->si.resize(numX * numY);
            this->li.resize(numX * numY);
            this->precon.resize(numX * numY);
            this->z.resize(numX * numY);
            this->search.resize(numX * numY);

            this->direction.resize(numX * numY);
            this->residual.resize(numX * numY);
            this->pressure.resize(numX * numY);
            this->residual.resize(numX * numY);
            this->Ad.resize(numX * numY);

            this->collisions.resize(numParticles);
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
            this->vaCopy.resize(numParticles * 4);
            texture.loadFromFile("white_circle.png");
            texture.generateMipmap();
            auto const texture_size = static_cast<sf::Vector2f>(texture.getSize());
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

            textureSizeX = texture_size.x;
            textureSizeY = texture_size.y;

            this->scalingFactor = 2 * radius;

            this->scaledWIDTH = std::ceil(static_cast<float>(WIDTH) / scalingFactor);
            this->scaledHEIGHT = std::ceil(static_cast<float>(HEIGHT) / scalingFactor);

            collisionGrid = CollisionGrid(scaledWIDTH, scaledHEIGHT);

            cellOccupantsGrid = CollisionGrid(numX, numY);

            // initializing particle positions
            float separation = radius / separationInit;  // gridsize: 65: 2.5; 70: 2.3; 90: 2.f; 130: 1.2

            int wideNum = std::floor((WIDTH - 2 * cellSpacing - 2) / (radius * separation));
            int highNum = numParticles / wideNum;

            float starting_px = radius + cellSpacing + 2;//(WIDTH - (radius * separation * wideNum)) / 2 + radius;
            float starting_py = (HEIGHT - (radius * separation * highNum)) / 2 + radius;

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

                px += this->radius * separation;

                if ((i + 1) % wideNum == 0) {
                    px = starting_px;
                    if (offset) {
                        px += this->radius;
                    }
                    py += this->radius * separation;
                    offset = !offset;
                }
            }

            n = this->numY;
            for (int i = 0; i < this->numX; ++i) {
                for (int j = 0; j < this->numY; ++j) {
                    if (i == 0 || j == 0 || i == this->numX - 1 || j == this->numY - 1) {
                        this->cellType[i * n + j] = SOLID_CELL;
                    }
                }
            }

            this->moveDist = 2 * radius;
            this->checkseparationDist = moveDist * moveDist;

            this->objectRadius = 50;//cellSpacing * 3;
            this->objectX = std::floor(WIDTH / 2);
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

            // linearly interpolate between the values in colorMap to create a gradient array 
            float num_colors = velColorMap.size() - 1; // number of colors - 1
            float num_steps = 1.f * velGradient.size() / num_colors; //num_steps = 50 * key_range
            int index = 0;
            for (int i = 0; i < num_colors; ++i) {  
                for (int x = 0; x < num_steps; ++x) {
                    float t = 1.f * x / num_steps;  // Interpolation factor
                    // Linear interpolation for r, g, b values between colorMap[i] andcolorMap [i+1]
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

            /*this->spacing = 2 * this->radius;
            this->tableSize = 2 * numParticles;
            this->cellCount.resize(this->tableSize + 1);
            this->particleArray.resize(numParticles);
            this->allHashCells.resize(numParticles);*/
    }

    void drawCells(sf::RenderWindow& window) {
        for (int i = 0; i < numX; ++i) {
            for (int j = 0; j < numY; ++j) {
                int idx = i * numY + j;
                if (cellType[idx] == SOLID_CELL) {
                    cellDrawer.setPosition(i * cellSpacing, j * cellSpacing);
                    cellDrawer.setFillColor(sf::Color(cellColor[3 * idx], cellColor[3 * idx + 1], cellColor[3 * idx + 2]));
                    window.draw(cellDrawer);
                }
            }
        }
    }

    float clamp(float x, float min, float max) {
        return (x < min) * min + (x > max) * max + (x >= min && x <= max) * x;
    }

    float sign(float x) {
        return (x < 0) * -1.f + (x > 0) * 1.f;
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

    void integrate(const uint32_t startIndex, const uint32_t endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            this->positions[2 * i] += this->velocities[2 * i] * dt;
            this->positions[2 * i + 1] += this->velocities[2 * i + 1] * dt;
            this->velocities[2 * i + 1] += gravityX * dt;
            this->velocities[2 * i] += gravityY * dt;

            if (this->fireActive && this->renderPattern == 3 && positions[2 * i + 1] < HEIGHT - cellSpacing - 10) {
                this->velocities[2 * i + 1] -= fireStrength * temperatures[i] * dt;
                if (this->temperatures[i] > 0) {
                    this->temperatures[i] -= tempDiffusion;
                }
                if (collisions[i] == 0) {
                    this->temperatures[i] = 0;
                }
            }
        }
    }

    void integrateMulti() {
        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->integrate(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
            });
        }

        this->integrate(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();
    }

    void addObjectsToGrids()
    {
        collisionGrid.clear();
        cellOccupantsGrid.clear();

        const float minX = cellSpacing;
        const float maxX = WIDTH - cellSpacing;
        const float minY = cellSpacing;
        const float maxY = HEIGHT - cellSpacing;

        uint32_t i{0};
        for (int32_t index = 0; index < numParticles; ++index) {
            float x = positions[2 * index];
            float y = positions[2 * index + 1];
            if (x > minX && x < maxX && y > minY && y < maxY) {

                int32_t gridX = x / scalingFactor;
                int32_t gridY = y / scalingFactor;
                collisionGrid.addAtom(gridX, gridY, i);
                
                // THIS SHIT IS TAKING A LONG FUCKING TIME
                int32_t cellOccupantsX = x / cellSpacing;
                int32_t cellOccupantsY = y / cellSpacing;
                cellOccupantsGrid.addAtom(cellOccupantsX, cellOccupantsY, i);

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

        if (dist2 < checkseparationDist && dist2 > eps) {
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

                    float dx = positions[particleIndex * 2] - mouseX;
                    float dy = positions[particleIndex * 2 + 1] - mouseY;
                    float d2 = dx * dx + dy * dy;

                    if (d2 > checkForceObjectseparationDist || d2 == 0.0f) continue;

                    float d = std::sqrt(d2);

                    float edgeT = d / forceObjectRadius;
            
                    float centerT = 1 - edgeT;

                    velocities[2 * particleIndex] += (dx * strength - velocities[2 * particleIndex]) * centerT * dt;
                    velocities[2 * particleIndex + 1] += (dy * strength - velocities[2 * particleIndex + 1]) * centerT * dt;
                }
            }
        }
    }

    void solveCollisions()
    {
        const uint32_t thread_count = thread_pool.m_thread_count;
        const uint32_t slice_count  = thread_count * 2;
        const uint32_t slice_size   = (collisionGrid.width / slice_count) * collisionGrid.height;
        const uint32_t last_cell    = 2 * thread_count * slice_size;
        
        for (uint32_t i{0}; i < thread_count; ++i) {
            thread_pool.addTask([this, i, slice_size]{
                uint32_t const start{2 * i * slice_size};
                uint32_t const end  {start + slice_size};
                solveCollisionThreaded(start, end);
            });
        }
        
        if (last_cell < collisionGrid.data.size()) {
            thread_pool.addTask([this, last_cell]{
                solveCollisionThreaded(last_cell, static_cast<uint32_t>(collisionGrid.data.size()));
            });
        }
        thread_pool.waitForCompletion();
        
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
                const float remove = 10.f; // 5
                if (this->renderPattern == 3 && this->temperatures[i] < tempgradient.size() && this->positions[2 * i] > WIDTH / remove && this->positions[2 * i] < WIDTH - WIDTH / remove) {
                    this->temperatures[i] += groundConductivity;
                }
            }
        }
    }

    void constrainWallsMulti() {
        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->constrainWalls(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
            });
        }

        this->constrainWalls(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();
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
    
        if (pdx > 0 && pdy > 0 && cellType[idx + dirX] != SOLID_CELL && cellType[idx + dirY] != SOLID_CELL) {
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
        int localX = px / cellSpacing;
        int localY = py / cellSpacing;

        int x0 = std::max(0, localX - 1);
        int x1 = std::min(numX - 1, localX + 1);
        int y0 = std::max(0, localY - 1);
        int y1 = std::min(numY - 1, localY + 1);

        float prevX = px;
        float prevY = py;

        for (int i = x0; i <= x1; ++i) {
            for (int j = y0; j <= y1; ++j) {
                if (cellType[i * n + j] == SOLID_CELL) {
                    float nx = 0;
                    float ny = 0;
                    float dist = calculateBoxNormals(i * cellSpacing + halfSpacing, j * cellSpacing + halfSpacing, px, py, nx, ny);
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
            separateParticle(positions[2 * i], positions[2 * i + 1], velocities[2 * i], velocities[2 * i + 1]);
        }
    }

    void collideSurfacesMulti() {
        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->collideSurfaces(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
            });
        }

        this->collideSurfaces(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();
    }

    void drawSolids() {
        int localX = static_cast<int>(mouseX / cellSpacing);
        int localY = static_cast<int>(mouseY / cellSpacing);

        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                int idx = (localX + i) * numY + localY + j;
                if (cellType[idx] != SOLID_CELL && localX + i % numX > 0 && localY + j > 0 && localX + i % numX < numX - 1 && localY + j < numY - 1) {
                    cellType[idx] = SOLID_CELL;
                }
            }
        }
    }

    void eraseSolids() {
        int localX = static_cast<int>(mouseX / cellSpacing);
        int localY = static_cast<int>(mouseY / cellSpacing);

        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                int idx = (localX + i) * numY + localY + j;
                if (cellType[idx] == SOLID_CELL && localX + i % numX > 0 && localY + j > 0 && localX + i % numX < numX - 1 && localY + j < numY - 1) {
                    cellType[idx] = AIR_CELL;
                }
            }
        }
    }

    void cacheTransferNodes(int32_t start, int32_t end, float halfHeight, int32_t component) {
        const float h2 = halfHeight;

        const float dx = (component != 0) * h2;
        const float dy = (component == 0) * h2;

        for (int32_t i = start; i < end; ++i) {
            float x = this->positions[2 * i];
            float y = this->positions[2 * i + 1];
            x = this->clamp(x, cellSpacing, (this->numX - 1) * cellSpacing);
            y = this->clamp(y, cellSpacing, (this->numY - 1) * cellSpacing);
            // x0 is the grid position to the left of the particle, x1 is the position to the right of the particle. Both can only go up to the second to last cell to the right in the grid because we dont want to be changing wall velocities
            int x0 = std::max(1, std::min(static_cast<int>(std::floor((x - dx) * invSpacing)), this->numX - 2)); // - 1
            // basically x - xCell to get the weight of that cell in relation to the particle
            // in this case, x is moved over to x - dx, and xCell is just grid position of x multiplied by grid spacing
            float tx = ((x - dx) - x0 * cellSpacing) * invSpacing;
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
            float ty = ((y - dy) - y0 * cellSpacing) * invSpacing;
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

        const int32_t numThreads = thread_pool.m_thread_count;
        const int32_t numParticlesPerThread = numParticles / numThreads;
        const int32_t numMissedParticles = numParticles - numThreads * numParticlesPerThread;

        for (int32_t i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i](){
                cacheTransferNodes(numParticlesPerThread * i, numParticlesPerThread * i + numParticlesPerThread, halfSpacing, 0);
                cacheTransferNodes(numParticlesPerThread * i, numParticlesPerThread * i + numParticlesPerThread, halfSpacing, 1);
            });
        }

        cacheTransferNodes(numParticles - numMissedParticles, numParticles, halfSpacing, 0);
        cacheTransferNodes(numParticles - numMissedParticles, numParticles, halfSpacing, 1);

        thread_pool.waitForCompletion();
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
        const float h2 = 0.5 * cellSpacing;

        std::copy(std::begin(this->u), std::end(this->u), std::begin(this->prevU));
        std::copy(std::begin(this->v), std::end(this->v), std::begin(this->prevV));
        std::fill(begin(this->du), end(this->du), 0.f);
        std::fill(begin(this->dv), end(this->dv), 0.f);
        std::fill(begin(this->u), end(this->u), 0.f);
        std::fill(begin(this->v), end(this->v), 0.f);

        // initialize outside cells to solid, every inside cell to air
        for (int i = 0; i < this->numX; ++i) {
            for (int j = 0; j < this->numY; ++j) {
                if (this->cellType[i * n + j] == SOLID_CELL) {
                    this->cellColor[3 * (i * n + j)] = 100;
                    this->cellColor[3 * (i * n + j) + 1] = 100;
                    this->cellColor[3 * (i * n + j) + 2] = 100;
                }
                else if (this->cellType[i * n + j] != SOLID_CELL) {
                    this->cellType[i * n + j] = AIR_CELL;
                       
                    /*this->cellColor[3 * (i * n + j)] = 0;
                    this->cellColor[3 * (i * n + j) + 1] = 250;
                    this->cellColor[3 * (i * n + j) + 2] = 255;*/
                        
                }
            }
        }
        // initialize all cells that particles are in to fluid
        for (int i = 0; i < this->numParticles; ++i) {
            float x = this->positions[2 * i];
            float y = this->positions[2 * i + 1];

            int xi = this->clamp(std::floor(x * invSpacing), 0, this->numX - 1);
            int yi = this->clamp(std::floor(y * invSpacing), 0, this->numY - 1);

            int cellNr = xi * n + yi;
            if (this->cellType[cellNr] == AIR_CELL) {
                this->cellType[cellNr] = FLUID_CELL;
        
                /*this->cellColor[3 * cellNr] = 0;
                this->cellColor[3 * cellNr + 1] = 150;
                this->cellColor[3 * cellNr + 2] = 255;*/
                
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

            const float pvx = this->velocities[ui];
            const float pvy = this->velocities[vi];
           
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
                this->velocities[ui] = (1.f - flipRatio) * picV + flipRatio * flipV;
            }

            if (divY > 0.f) {
                picV = (valid0v * d0_v * this->v[nr0_v] + valid1v * d1_v * this->v[nr1_v] + valid2v * d2_v * this->v[nr2_v] + valid3v * d3_v * this->v[nr3_v]) / divY;

                corr = (valid0v * d0_v * (this->v[nr0_v] - this->prevV[nr0_v]) + valid1v * d1_v * (this->v[nr1_v] - this->prevV[nr1_v]) + valid2v * d2_v * (this->v[nr2_v] - this->prevV[nr2_v]) + valid3v * d3_v * (this->v[nr3_v] - this->prevV[nr3_v])) / divY;
                flipV = pvy + corr;
                this->velocities[vi] = (1.f - flipRatio) * picV + flipRatio * flipV;
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

            const float pvx = this->velocities[ui];
            const float pvy = this->velocities[vi];
           
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

            const float pvx = this->velocities[ui];
            const float pvy = this->velocities[vi];
    
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
            thread_pool.addTask([&, start, end]() {
                this->transferToUGrid(start, end);
            });
            thread_pool.addTask([&, start, end]() {
                this->transferToVGrid(start, end);
            });
        }

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

        std::fill(begin(this->particleDensity), end(this->particleDensity), 0.f);

        for (int i = 0; i < numParticles; ++i) {
            float x = this->positions[2 * i];
            float y = this->positions[2 * i + 1];

            x = this->clamp(x, cellSpacing, (this->numX - 1) * cellSpacing);
            y = this->clamp(y, cellSpacing, (this->numY - 1) * cellSpacing);

            int x0 = std::max(1, std::min((int)(std::floor((x - halfSpacing) * invSpacing)), this->numX - 2));
            float tx = ((x - halfSpacing) - x0 * cellSpacing) * invSpacing;
            int x1 = std::min(x0 + 1, this->numX - 1);

            int y0 = std::max(1, std::min((int)(std::floor((y - halfSpacing) * invSpacing)), this->numY - 2));
            float ty = ((y - halfSpacing) - y0 * cellSpacing) * invSpacing;
            int y1 = std::min(y0 + 1, this->numY - 2);
           
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

    void setUpResidual() {
        double scale = 1.0 / cellSpacing;

        std::fill(begin(residual), end(residual), 0.0);

        for (int32_t i = 1; i < numX - 1; ++i) {
            for (int32_t j = 1; j < numY - 1; ++j) {
                int32_t idx = i * n + j;
                if (cellType[idx] != FLUID_CELL)  {
                    continue;
                }

                float divergence = -(this->u[idx + n] - this->u[idx] + this->v[idx + 1] - this->v[idx]);

                if (this->particleRestDensity > 0.f) {
                    float compression = this->particleDensity[idx] - this->particleRestDensity;
                    divergence += (compression > 0.f) * this->k * compression;
                }

                residual[idx] = divergence;

            }
        }
    }

    void ScaledAdd(std::vector<double>& a, std::vector<double>& b, double c) {
        for (int i = 0; i < numX * numY; ++i) {
            if (cellType[i] == FLUID_CELL) {
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
            thread_pool.addTask([&, i, start, end]() {
                this->Dot(a, b, start, end, dotProducts[i]);
            });
        }

        thread_pool.waitForCompletion();

        double res = 0.0;
        for (double el : dotProducts) {
            res += el;
        }

        return res;
    }

    void Dot(std::vector<double> *a, std::vector<double> *b, int start, int end, double& res) {
        for (int i = start; i < end; ++i) {
            if (cellType[i] == FLUID_CELL) {
                res += (*a)[i] * (*b)[i];
            }
        }
    }

    void EqualsPlusTimesMulti(std::vector<double> *a, std::vector<double> *b, double c) {
        const int32_t numThreads_ = 1;
        const int32_t numColumnsPerThread = (numX * numY) / numThreads_;
        const int32_t numMissedColumns = numX * numY - numColumnsPerThread * numThreads_;

        for (int i = 0; i < numThreads_; ++i) {
            thread_pool.addTask([&, i]() {
                this->EqualsPlusTimes(a, b, c, i * numColumnsPerThread, i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->EqualsPlusTimes(a, b, c, numX * numY - numMissedColumns, numX * numY);

        thread_pool.waitForCompletion();
    }

    void EqualsPlusTimes(std::vector<double> *a, std::vector<double> *b, double c, int start, int end) {
        for (int i = start; i < end; ++i) {
            if (cellType[i] == FLUID_CELL) {
                (*a)[i] = (*b)[i] + (*a)[i] * c;
            }
        }
    }

    void applyPressureMulti() {
        const int32_t numThreads_ = 1;
        const int32_t numColumnsPerThread = (numX - 2) / numThreads_;
        const int32_t numMissedColumns = (numX - 2) - numColumnsPerThread * numThreads_;

        for (int i = 0; i < numThreads_; ++i) {
            thread_pool.addTask([&, i]() {
                this->applyPressure(1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->applyPressure(numX - 1 - numMissedColumns, numX - 1);

        thread_pool.waitForCompletion();
    }

    void applyPressure(int start, int end) {
        //const float density = 1000.f;
        //const float scale = dt / (density * cellSpacing);
    
        for (int i = start; i < end; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int idx = i * n + j;
                int leftIdx =  idx - n;
                int upIdx = idx - 1;

                if ((cellType[idx] == FLUID_CELL || cellType[leftIdx] == FLUID_CELL) &&
                    cellType[idx] != SOLID_CELL && cellType[leftIdx] != SOLID_CELL) {
                    float p = pressure[idx];
                    float pLeft = pressure[leftIdx];
                    u[idx] -= 1 * (p - pLeft);
                }
    
                if ((cellType[idx] == FLUID_CELL || cellType[upIdx] == FLUID_CELL) &&
                    cellType[idx] != SOLID_CELL && cellType[upIdx] != SOLID_CELL) {
                    float p = pressure[idx];
                    float pTop = pressure[upIdx];
                    v[idx] -= 1 * (p - pTop);
                }
            }
        }
    }

    void setUpA() {
        float scale = 1;//dt / (cellSpacing * cellSpacing); // also see what happens when you divide by 1000 (density of water) as well

        std::fill(Adiag.begin(), Adiag.end(), 0.0);

        // Ax[i] is the cell to the right, Ay[i] is the cell below
        std::fill(si.begin(), si.end(), 0.0);
        std::fill(li.begin(), li.end(), 0.0);

        for (int i = 1; i < numX - 1; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int idx = i * n + j;

                if (cellType[idx] != FLUID_CELL) continue;

                int left = cellType[idx - n];
                int right = cellType[idx + n];
                int bottom = cellType[idx + 1];
                int top = cellType[idx - 1];

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

                if (cellType[idx] != FLUID_CELL) continue;

                double e = Adiag[idx];

                // si is short for swipe iteration, just some slang for you new gens
                    // si[index] = 0 if the cell in the next swipe iteration touching the current is not fluid, else si[index] = -1

                // li is short for level iteration, try and keep up with the slang
                    // li[index] = 0 if the cell in the next level iteration touching the current cell is not fluid, else li[index] = -1

                // if I add/subtract something by S in a comment, that means I'm taking the value of it at the next/previous swipe iteration
                // same for level iteration, denoted with L

                // if cell - S is fluid, then do: (si - S) * (precon - S), and (li - S) * (precon - S)
                if (cellType[idx - 1] == FLUID_CELL) {
                    double px = si[idx - 1] * precon[idx - 1];
                    double py = li[idx - 1] * precon[idx - 1];
                    e -= (px * px + tau * px * py);
                }

                // if cell - L is fluid, then do: (si - L) * (precon - L), and (li - L) * (precon - L)
                if (cellType[idx - n] == FLUID_CELL) {
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
                if (cellType[idx] != FLUID_CELL) continue;

                double t = (*a)[idx];

                // if cell - S is fluid, then do: (si - S) * (precon - S) * (dst - S)
                if (cellType[idx - 1] == FLUID_CELL) {
                    t -= si[idx - 1] * precon[idx - 1] * (*dst)[idx - 1];
                }

                // if cell - L is fluid, then do: (si - L) * (precon - L) * (dst - L)
                if (cellType[idx - n] == FLUID_CELL) {
                    t -= li[idx - n] * precon[idx - n] * (*dst)[idx - n];
                }

                (*dst)[idx] = t * precon[idx];
            }
        }

        for (int i = numX - 2; i > 0; --i) {
            for (int j = numY - 2; j > 0; --j) {
                int32_t idx = i * n + j;
                if (cellType[idx] != FLUID_CELL) continue;

                double t = (*dst)[idx];

                // if cell + S is fluid, then do: (si) * (precon) * (dst + S)
                if (cellType[idx + 1] == FLUID_CELL) {
                    t -= si[idx] * precon[idx] * (*dst)[idx + 1];
                }

                // if cell + L is fluid, then do: (li) * (precon) * (dst + L)
                if (cellType[idx + n] == FLUID_CELL) {
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

        std::fill(begin(this->p), end(this->p), 0.f);

        std::copy(std::begin(this->u), std::end(this->u), std::begin(this->prevU));
        std::copy(std::begin(this->v), std::end(this->v), std::begin(this->prevV));

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
            applyPreconditioner(&z, &residual);
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
        std::fill(begin(this->p), end(this->p), 0.f);

        std::copy(std::begin(this->u), std::end(this->u), std::begin(this->prevU));
        std::copy(std::begin(this->v), std::end(this->v), std::begin(this->prevV));

        std::fill(begin(pressure), end(pressure), 0.f);

        for (int iter = 0; iter < numPressureIters; ++iter) {
            for (int i = 1; i < this->numX - 1; ++i) {
                for (int j = 1; j < this->numY - 1; ++j) {
                    if (this->cellType[i * n + j] != FLUID_CELL) continue;

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
        }
    }

    float curl(int i, int j) {
        const int32_t n = this->numY;
        const float denom = 1.f / (2.f * cellSpacing);
        const int32_t leftType = cellType[(i - 1) * n + j] == FLUID_CELL;
        const int32_t rightType = cellType[(i + 1) * n + j] == FLUID_CELL;
        const int32_t topType = cellType[i * n + j - 1] == FLUID_CELL;
        const int32_t bottomType = cellType[i * n + j + 1] == FLUID_CELL;
        if (!leftType || !rightType || !topType || !bottomType) {
            return 0.f;
        }
        return ((this->v[(i + 1) * n + j] * bottomType - this->v[(i - 1) * n + j] * topType) - (this->u[i * n + j + 1] * rightType - this->u[i * n + j - 1] * leftType)) * denom;
    }

    float calcVorticitySquared(int i, int j) {
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
        const float extend = 20 * cellSpacing;
        if (mouseDown) {
            float vx = (objectX - objectPrevX) * 100;
            float vy = (objectY - objectPrevY) * 100;
            for (int i = 1; i < numX - 1; i++) {
                for (int j = 1; j < numY - 1; j++) {
                    int cellNr = i * n + j;
                    if (cellType[i * n + j] == SOLID_CELL) continue;
                    float dx = (i + 0.5) * cellSpacing - objectX;
                    float dy = (j + 0.5) * cellSpacing - objectY;

                    if (dx * dx + dy * dy < objectRadius * objectRadius + extend) {
                        cellType[i * n + j] = FLUID_CELL;

                        this->cellColor[3 * cellNr] = 0;
                        this->cellColor[3 * cellNr + 1] = 150;
                        this->cellColor[3 * cellNr + 2] = 255;

                        this->u[cellNr] = vx;
                        this->u[cellNr + n] = vx;
                        this->v[cellNr] = vy;
                        this->v[cellNr + 1] = vy;
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

    void generate() {
        float separation = radius / separationInit;
        int wideNum = std::floor((2 * generatorRadius) / (radius * separation));
        int highNum = wideNum;
        int numAdded = wideNum * highNum;

        float starting_px = mouseX - generatorRadius + radius;
        float starting_py = mouseY - generatorRadius + radius;
        float px = starting_px;
        float py = starting_py;
        bool offset = true;

        std::vector<float> addToPositions;

        for (int i = 0; i < numAdded; ++i) {
            int cellX = px / cellSpacing;
            int cellY = py / cellSpacing;
            int cellNr = cellX * n + cellY;
            float prevPx = px;
            float prevPy = py;

            px += this->radius * separation;
            if ((i + 1) % wideNum == 0) {
                px = starting_px;
                if (offset) {
                    px += this->radius * separation;
                }
                py += this->radius * separation;
                offset = !offset;
            }

            if (cellType[cellNr] != AIR_CELL || px - radius < 0 || px + radius > WIDTH || py - radius < 0 || py + radius > HEIGHT) {
                continue;
            }

            addToPositions.push_back(prevPx);
            addToPositions.push_back(prevPy);
        }

        int addedTo = addToPositions.size() / 2;

        numParticles += addedTo;

        this->positions.resize(2 * numParticles);
        this->velocities.resize(2 * numParticles);
        this->particleColors.resize(3 * numParticles);
        this->va.resize(4 * numParticles);

        int start = numParticles - addedTo;

        for (int i = start; i < numParticles; i++) {
            int idx1 = 2 * i;
            int idx2 = 3 * i;
            int idx3 = 4 * i;
            int idx4 = i - start;

            positions[idx1] = addToPositions[idx4 * 2];
            positions[idx1 + 1] = addToPositions[idx4 * 2 + 1];
            
            velocities[idx1] = 0.f;
            velocities[idx1 + 1] = 0.f;

            particleColors[idx2] = 255;
            particleColors[idx2 + 1] = 255;
            particleColors[idx2 + 2] = 255;

            va[idx3].texCoords = {0.f, 0.f};
            va[idx3 + 1].texCoords = {textureSizeX, 0.f};
            va[idx3 + 2].texCoords = {textureSizeX, textureSizeY};
            va[idx3 + 3].texCoords = {0.f, textureSizeY};
        }

        particlesPerThread = numParticles / numThreads;
        numMissedParticles = numParticles - numThreads * particlesPerThread;

        this->collisions.resize(numParticles);
        this->temperatures.resize(numParticles);

        this->nr0.resize(2 * numParticles);
        this->nr1.resize(2 * numParticles);
        this->nr2.resize(2 * numParticles);
        this->nr3.resize(2 * numParticles);

        this->d0.resize(2 * numParticles);
        this->d1.resize(2 * numParticles);
        this->d2.resize(2 * numParticles);
        this->d3.resize(2 * numParticles);
    }

    void remove() {
        const int32_t numCovered = std::ceil(generatorRadius / scalingFactor);
        const uint32_t mouseColumn = std::floor(mouseX / scalingFactor);
        const uint32_t mouseRow = std::floor(mouseY / scalingFactor);
        const size_t double_len = positions.size();
        const size_t triple_len = particleColors.size();
        const size_t quadruple_len = 2 * double_len;
    
        size_t vaSize = va.getVertexCount();
        vaCopy.resize(vaSize);
    
        for (int i = 0; i < vaSize; ++i) {
            vaCopy[i] = va[i];
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
                        positions[doubleid] = positions[double_len - 2 - doubleremove];
                        positions[doubleid + 1] = positions[double_len - 1 - doubleremove];

                        velocities[doubleid] = velocities[double_len - 2 - doubleremove];
                        velocities[doubleid + 1] = velocities[double_len - 1 - doubleremove];

                        particleColors[tripleid] = particleColors[triple_len - 3 - tripleremove];
                        particleColors[tripleid + 1] = particleColors[triple_len - 2 - tripleremove];
                        particleColors[tripleid + 2] = particleColors[triple_len - 1 - tripleremove];

                        vaCopy[quadrupleid] = vaCopy[quadruple_len - 4 - quadrupleremove];
                        vaCopy[quadrupleid + 1] = vaCopy[quadruple_len - 3 - quadrupleremove];
                        vaCopy[quadrupleid + 2] = vaCopy[quadruple_len - 2 - quadrupleremove];
                        vaCopy[quadrupleid + 3] = vaCopy[quadruple_len - 1 - quadrupleremove];
                    }

                    remove++;
                }
            }
        }

        numParticles -= remove;
    
        this->positions.resize(2 * numParticles);
        this->velocities.resize(2 * numParticles);
        this->particleColors.resize(3 * numParticles);
        this->vaCopy.resize(4 * numParticles);
    
        va.resize(4 * numParticles);
        for (int i = 0; i < va.getVertexCount(); ++i) {
            va[i] = vaCopy[i];
        }
    
        particlesPerThread = numParticles / numThreads;
        numMissedParticles = numParticles - numThreads * particlesPerThread;
    
        this->collisions.resize(numParticles);
        this->temperatures.resize(numParticles);
    
        this->nr0.resize(2 * numParticles);
        this->nr1.resize(2 * numParticles);
        this->nr2.resize(2 * numParticles);
        this->nr3.resize(2 * numParticles);
    
        this->d0.resize(2 * numParticles);
        this->d1.resize(2 * numParticles);
        this->d2.resize(2 * numParticles);
        this->d3.resize(2 * numParticles);
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
        int localX = static_cast<int>(mouseX / cellSpacing);
        int localY = static_cast<int>(mouseY / cellSpacing);
        int idx = localX * n + localY;

        float drawPosX = localX * cellSpacing + halfSpacing;
        float drawPosY = localY * cellSpacing + halfSpacing;

        if (leftMouseDown || !rightMouseDown) {
            pencil.setFillColor(sf::Color(0, 150, 0));
        }
        else {
            pencil.setFillColor(sf::Color(150, 0, 0));
        }

        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                if (localX + i > 0 && localY + j > 0 && localX + i < numX - 1 && localY + j < numY - 1) {
                    pencil.setPosition(drawPosX + i * cellSpacing, drawPosY + j * cellSpacing);
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
                if (cellType[idx] == FLUID_CELL) {
                    float div = fabsf(u[(i + 1) * n + j] - u[idx] + v[idx + 1] - v[idx]);

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
                if (cellType[idx] != FLUID_CELL) continue;
                float div = fabsf(u[(i + 1) * n + j] - u[idx] + v[idx + 1] - v[idx]);
                cellDrawer.setPosition(i * cellSpacing, j * cellSpacing);

                div /= maxDiv;

                int redScale = 255 * (div);
                int greenScale = 255 * (1.f - div);

                cellDrawer.setFillColor(sf::Color(redScale, greenScale, 0));
                window.draw(cellDrawer);
            }
        }
    }

    void updateVertexArrayVelocity(uint32_t startIndex, uint32_t endIndex) {
        int32_t velGradientSize = velGradient.size() - 1;
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            int i = 4 * index;
            const float px = positions[2 * index];
            const float py = positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            sf::Color color;

            int vel = (int)(velocities[2 * index] * velocities[2 * index] + velocities[2 * index + 1] * velocities[2 * index + 1]) / 15000; 
            
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
            const float px = positions[2 * index];
            const float py = positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            int cellOccupantsX = px / cellSpacing;
            int cellOccupantsY = py / cellSpacing;

            float vort = static_cast<float>(std::min(255, static_cast<int>(calcVorticitySquared(cellOccupantsX, cellOccupantsY) * 5))); 
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
            const float s = 0.01f;

            int i = 4 * index;
            const float px = positions[2 * index];
            const float py = positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            particleColors[3 * index] = clamp(this->particleColors[3 * index] - s, 0, 255);
            particleColors[3 * index + 1] = clamp(this->particleColors  [3 * index + 1] - s, 0, 255);
            particleColors[3 * index + 2] = clamp(this->particleColors  [3 * index + 2] + s, 0, 255);

            const int xi = clamp(std::floor(px * invSpacing), 1,    this->numX - 1);
            const int yi = clamp(std::floor(py * invSpacing), 1,    this->numY - 1);
            const int cellNr = xi * this->numY + yi;

            const float d0 = this->particleRestDensity;

            if (d0 > 0.f) {
                const float relDensity = this->particleDensity  [cellNr] / d0;
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

    void drawVelocityMulti() {
        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->updateVertexArrayVelocity(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
            });
        }

        this->updateVertexArrayVelocity(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();
    }

    void drawVorticityMulti() {
        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->updateVertexArrayVorticity(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
            });
        }

        this->updateVertexArrayVorticity(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();
    }

    void drawDiffusionMulti() {
        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->updateVertexArrayDiffusion(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
            });
        }

        this->updateVertexArrayDiffusion(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();
    }

    void drawTemperatureMulti() {
        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->updateVertexArrayTemperature(i * particlesPerThread, i * particlesPerThread + particlesPerThread);
            });
        }

        this->updateVertexArrayTemperature(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();
    }

    void drawParticlesVertex(sf::RenderWindow& window) {
        window.draw(va, states);
    }

    void addValueToAverage(float& value, float newValue) {
        value += (newValue - value) / steps;
    }

    void update(float dt_, sf::RenderWindow& window, bool leftMouseDown, bool justPressed, bool rightMouseDown) {
        auto start = std::chrono::high_resolution_clock::now();
        sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
        this->mouseX = mouse_pos.x;
        this->mouseY = mouse_pos.y;
        if (!stop || step) {
            this->simulate(dt_, leftMouseDown, justPressed, rightMouseDown);
            step = false;
        }

        this->render(window);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        addValueToAverage(SimStepTime, duration.count());
    }

    void simulate(float dt_, bool leftMouseDown_, bool rightMouseDown_, bool justPressed) {

        ++steps;

        // order of need of implementation/optimization:
            // 1) incompressibility -- optimize PCG
            // 2) implement stable volume correction
            // 3) updateDensity -- Same idea as to grid

        // collision, rendering, and to particles all great
        auto start = std::chrono::high_resolution_clock::now();

        dt = dt_;

        this->integrateMulti();

        if (fireActive) {
            std::fill(begin(collisions), end(collisions), 0);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count());





        start = std::chrono::high_resolution_clock::now();

        addObjectsToGrids();

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(FillGridTime, duration.count());





        start = std::chrono::high_resolution_clock::now();
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

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count());





        start = std::chrono::high_resolution_clock::now();
        solveCollisions();
        
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(CollisionTime, duration.count());





        start = std::chrono::high_resolution_clock::now();

        this->collideSurfacesMulti();

        this->constrainWallsMulti();

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ObstacleCollisionTime, duration.count());





        start = std::chrono::high_resolution_clock::now();
        this->transfer(true);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ToGridTime, duration.count());



        start = std::chrono::high_resolution_clock::now();
        this->updateParticleDensity();
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(DensityUpdateTime, duration.count());



        start = std::chrono::high_resolution_clock::now();
        if (rigidObjectActive) {
            this->includeRigidObject(leftMouseDown, justPressed);
        }


        if (vorticityStrength > 0) {
            this->applyVorticityConfinementRedBlack();
        }
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count());

        start = std::chrono::high_resolution_clock::now();
        this->PCGproject();
        //this->SORproject();
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ProjectionTime, duration.count());




        start = std::chrono::high_resolution_clock::now();
        this->transfer(false);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ToParticlesTime, duration.count());
    }

    void render(sf::RenderWindow& window) {
        auto start = std::chrono::high_resolution_clock::now();
        this->drawCells(window);

        if (renderPattern == 0) {
            this->drawDiffusionMulti();
        }

        else if (renderPattern == 1) {
            this->drawVelocityMulti();
        }

        else if (renderPattern == 2) {
            this->drawVorticityMulti();
        }

        else if (renderPattern == 3) {
            this->drawTemperatureMulti();
        }

        else if (renderPattern == 4) {
            this->DrawDivergences(window);
        }

        if (renderPattern != 4) {
            this->drawParticlesVertex(window);
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
        this->gravityX += add;
    }

    float getGravityX() {
        return this->gravityX;
    }

    void addToGravityY(float add) {
        this->gravityY += add;
    }

    float getGravityY() {
        return this->gravityY;
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

    float getTimeForseparation() {
        return this->timeForseparation;
    }

    float getTimeForInc() {
        return this->timeForIncompressibility;
    }

    float getTimeForTrans() {
        return this->timeForTransfer;
    }

    void setNextRenderPattern() {
        this->renderPattern++;
        if (this->renderPattern > 4) {
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

    bool getStop() {
        return this->stop;
    }

    void setStop(bool set) {
        this->stop = set;
    }

    bool getStep() {
        return this->step;
    }

    void setStep(bool set) {
        this->step = set;
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
                float uX = cellSpacing + i * cellSpacing;
                float uY = 1.5 * cellSpacing + j * cellSpacing;
                line[0].position = sf::Vector2f(uX, uY);
                line[0].color  = sf::Color(255, 150, 0);
                line[1].position = sf::Vector2f(uX + u[i * n + j], uY);
                line[1].color = sf::Color(255, 150, 0);
                window.draw(line);

                //draw v lines (top bottom)
                float vX = 1.5 * cellSpacing + i * cellSpacing;
                float vY = cellSpacing + j * cellSpacing;
                line[0].position = sf::Vector2f(vX, vY);
                line[0].color  = sf::Color::Red;
                line[1].position = sf::Vector2f(vX, vY + v[i * n + j]);
                line[1].color = sf::Color::Red;
                window.draw(line);
            }
        }
    }
};
