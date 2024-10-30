#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <array>
#include <sstream>
#include <iomanip>

class Fluid {
    float density;
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
    int tableSize;
    int numObjects;
    std::vector<int> cellCount;
    std::vector<int> particleArray;
    std::vector<int> allHashCells;
    float spacing;
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
    sf::RectangleShape cellDrawer;
    sf::CircleShape uDrawer;
    sf::CircleShape vDrawer;
    std::vector<float> uColor;
    std::vector<float> vColor;

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

    std::vector<int> particleColors;

    float gravity;

    std::array<std::array<int, 3>, 100> gradient;
    std::array<std::array<int, 3>, 4> colorMap{{{0, 51, 102}, {0, 153, 204}, {102, 255, 204}, {255, 255, 255}}};
    // some nice gradients to put into colorMap:
    // scientific: {0, 150, 255}, {0, 255, 0}, {255, 255, 0}, {255, 0, 0}
    // night ocean: {0, 51, 102}, {0, 153, 204}, {102, 255, 204}, {255, 255, 255}
    // sunset: {0, 0, 64}, {128, 0, 128}, {255, 128, 0}, {255, 255, 0}
    // orange to white: {102, 51, 0}, {204, 122, 0}, {255, 153, 51}, {255, 255, 255}
    // ice: {0, 102, 204}, {173, 216, 230}, {224, 255, 255}, {255, 250, 250}
    // lava: {128, 0, 0}, {255, 69, 0}, {255, 140, 0}, {255, 215, 0}
    // deep space: {0, 0, 32}, {64, 0, 128}, {128, 0, 255}, {192, 192, 255}
    // dark blue: {0, 0, 128}, {0, 128, 255}, {255, 128, 0}, {255, 255, 0}
    // lightning mcqueen: {255, 0, 0}, {255, 69, 0}, {255, 165, 0}, {255, 255, 0}
    // rainbow: {255, 0, 0}, {255, 255, 0}, {0, 255, 0}, {0, 200, 255} 

    float nX;
    float nY;
    std::vector<int> cellCount2;
    float tableSize2;

    sf::VertexArray va{sf::PrimitiveType::Quads};

    sf::Transform tf;

    sf::Texture texture;

    sf::RenderStates states;

public:
    Fluid(float density, float WIDTH, float HEIGHT, float cellSpacing, int numParticles, float gravity)
        : numX(std::floor(WIDTH / cellSpacing)), numY(std::floor(HEIGHT / cellSpacing)), numCells(numX * numY), numParticles(numParticles), WIDTH(WIDTH), HEIGHT(HEIGHT), gravity(gravity) {
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
            auto const texture_size = static_cast<sf::Vector2f>(texture.getSize());
            for (int index = 0; index < numParticles; ++index) {
                int i = 4 * index;
                va[i].texCoords = {0.f, 0.f};
                va[i + 1].texCoords = {texture_size.x, 0.f};
                va[i + 2].texCoords = {texture_size.x, texture_size.y};
                va[i + 3].texCoords = {0.f, texture_size.y};
            }
            states.texture = &texture;

            this->cellSpacing = std::max(WIDTH / numX, HEIGHT / numY);
            this->invSpacing = 1.f / this->cellSpacing;

            this->radius = 0.3 * cellSpacing;
            this->particleDrawer.setOutlineThickness(0.f);
            this->particleDrawer.setRadius(this->radius);
            this->particleDrawer.setOrigin(this->radius, this->radius);
            this->particleDrawer.setFillColor(sf::Color(0, 150, 255));

            this->spacing = 2 * this->radius;
            this->tableSize = 2 * numParticles;
            this->cellCount.resize(this->tableSize + 1);
            this->particleArray.resize(numParticles);
            this->allHashCells.resize(numParticles);
            this->nY = std::ceil(1.f * HEIGHT / this->spacing);
            this->nX = std::ceil(1.f * WIDTH / this->spacing);
            this->tableSize2 = this->nX * this->nY;
            this->cellCount2.resize(tableSize2 + 1);


            // initialize particle positions
            int rowNum = std::floor(std::sqrt(numParticles));

            int seperation = 4;
            int starting_px = (WIDTH - (radius * seperation * rowNum)) / 2 + radius;
            int starting_py = (HEIGHT - (radius * seperation * rowNum)) / 2 + radius;

            int px = starting_px;
            int py = starting_py;

            int addTo = numParticles - rowNum * rowNum;

            bool offset = false;
            for (int i = 0; i < rowNum * rowNum + addTo; ++i) {
                this->positions[i * 2] = px;
                this->positions[i * 2 + 1] = py;
                this->particleColors[3 * i] = 0;
                this->particleColors[3 * i + 1] = 150;
                this->particleColors[3 * i + 2] = 255;

                px += this->radius * seperation;

                if ((i + 1) % rowNum == 0) {
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

            this->colors.resize(numParticles);
            std::fill(begin(this->colors), end(this->colors), sf::Color(0, 150, 255));
            this->moveDist = 2 * radius;
            this->checkSeperationDist = moveDist * moveDist;

            this->cellDrawer.setSize(sf::Vector2f(cellSpacing, cellSpacing));
            this->cellDrawer.setOutlineThickness(1.f);
            this->cellDrawer.setOutlineColor(sf::Color::Black);

            this->objectRadius = 50;
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

            // linearly interpolate between the values in colorMap to create a gradient array 
            float num_colors = colorMap.size() - 1; // number of colors - 1
            float num_steps = 1.f * gradient.size() / num_colors; //num_steps = 50 * key_range
            int index = 0;
            for (int i = 0; i < num_colors; ++i) {  // Iterate over adjacent color pairs
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
    }

private:

    int hashCoords(int xi, int yi) {
        int hash = (intCoord(xi) * 92837111) ^ (intCoord(yi) * 689287499);
        return std::abs(hash) % this->tableSize;
    }

    int intCoord(int coord) {
        // add 1 so that you dont get 0 for intCoord(xi) or intCoord(yi) in the hashCoords
        return std::floor(coord / this->spacing) + 1;
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

    int getCell(float x, float y) {
        int xi = std::floor(x / this->spacing);
        int yi = std::floor(y / this->spacing);
        return xi * this->nY + yi;
    }

    // Red-black gauss seidel implementation of solving incompressibility. 
    // May or may not converge faster than normal iteration, I just wanted to see if it would be viable for multithreading. 
    void solvePass(const float overRelaxation, const bool red) {
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
    }

public:

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

    void integrate(const float dt, const sf::RenderWindow& window) {
        for (int i = 0; i < this->numParticles; ++i) {
            //this->velocities[2 * i + 1] += gravity * dt;
            this->positions[2 * i] += this->velocities[2 * i] * dt;
            this->positions[2 * i + 1] += this->velocities[2 * i + 1] * dt;
            this->velocities[2 * i + 1] += gravity * dt;
        }

        sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
        this->mouseX = mouse_pos.x;
        this->mouseY = mouse_pos.y;
    }

    void initializeSHConstantMem() {

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

    void initializeSH() {
        std::fill(this->cellCount2.begin(), this->cellCount2.end(), 0);
        std::fill(this->particleArray.begin(), this->particleArray.end(), 0);
        // initialize cells in cellCount
        for (int i = 0; i < this->numParticles; ++i) {
            if (positions[2 * i] > 0 && positions[2 * i] < WIDTH && positions[2 * i + 1] > 0 && positions[2 * i + 1] < HEIGHT) {
                int cell = this->getCell(positions[2 * i], positions[2 * i + 1]);
                this->cellCount2[cell]++;
            }
        }
        
        // calc partial sum
        int sum = 0;
        for (int i = 0; i < tableSize2 + 1; ++i) {
            sum += this->cellCount2[i];
            this->cellCount2[i] = sum;
        }
        this->cellCount2[tableSize2] = sum;
       
        // fill particle array
        for (int i = 0; i < this->numParticles; ++i) {
            if (positions[2 * i] > 0 && positions[2 * i] < WIDTH && positions[2 * i + 1] > 0 && positions[2 * i + 1] < HEIGHT) {
                int cell = this->getCell(positions[2 * i], positions[2 * i + 1]);
                cellCount2[cell]--;
                particleArray[cellCount2[cell]] = i;
            }
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

    void makeParticleQueries(int index) {
        for (int i = -1; i < 2; ++i) {
            for (int j = -int(std::min(positions[2 * index + 1] / this->spacing, 1.f)); j < int(std::min(std::ceil(this->nY - positions[2 * index + 1] / this->spacing), 2.f)); ++j) {
                int cell = this->getCell(this->positions[index * 2] + i * this->spacing, this->positions[index * 2 + 1] + j * this->spacing);
                if (cell < 0 || cell > tableSize2 - 1) continue;
                int start = this->cellCount2[cell];
                int end = this->cellCount2[cell + 1];

                for (int p = start; p < end; ++p) {

                    int otherParticleID = this->particleArray[p];

                    if (otherParticleID == index) continue;

                    float dx = this->positions[otherParticleID * 2] - this->positions[index * 2];
                    float dy = this->positions[otherParticleID * 2 + 1] - this->positions[index * 2 + 1];
                    float d2 = dx * dx + dy * dy;
                    if (d2 > checkSeperationDist || d2 == 0.0) continue;
                    float d = std::sqrt(d2); 
                    float s = 0.5 * (moveDist - d) / d;
                    dx *= s;
                    dy *= s;
                    this->positions[2 * index] -= dx;
                    this->positions[2 * index + 1] -= dy;
                    this->positions[2 * otherParticleID] += dx;
                    this->positions[2 * otherParticleID + 1] += dy;
                }
            }
        }  
    }

    void makeParticleQueriesMulti(int startColumn, int endColumn) {
        for (int c = 1; c < nX - 1; ++c) {
            for (int r = 1; r < nY - 1; ++r) {
                
                // add cellSpacing to the row and column values because the sim is offset by that value
                int cell = r + c * nY;
                int firstStart = this->cellCount[cell];
                int firstEnd = this->cellCount[cell + 1];

                for (int particleKey = firstStart; particleKey < firstEnd; particleKey++) {
                    
                    /*int particleIndex = this->particleArray[particleKey];

                    particleColors[3 * particleIndex] = 0;
                    particleColors[3 * particleIndex + 1] = 255;
                    particleColors[3 * particleIndex + 2] = 0;

                    for (int i = -1; i < 2; ++i) {
                        for (int j = -int(std::min(positions[2 * particleIndex + 1] / this->spacing, 1.f)); j < int(std::min(std::ceil(this->nY - positions[2 * particleIndex + 1] / this->spacing), 2.f)); ++j) {
                            int otherCell = this->getCell(this->positions[particleIndex * 2] + i * this->spacing, this->positions[particleIndex * 2 + 1] + j * this->spacing);
                            if (otherCell < 0 || otherCell > tableSize2 - 1) continue;

                            int start = this->cellCount2[otherCell];
                            int end = this->cellCount2[otherCell + 1];

                            for (int otherParticleKey = start; otherParticleKey < end; ++otherParticleKey) {
                                // this line is executing
                                int otherParticleID = this->particleArray[otherParticleKey];

                                // stopping here every time
                                if (otherParticleID == particleIndex) continue;

                                float dx = this->positions[otherParticleID * 2] - this->positions[particleIndex * 2];
                                float dy = this->positions[otherParticleID * 2 + 1] - this->positions[particleIndex * 2 + 1];
                                float d2 = dx * dx + dy * dy;
                                if (d2 > checkSeperationDist || d2 == 0.0) continue;
                                float d = std::sqrt(d2); 
                                float s = 0.5 * (moveDist - d) / d;
                                dx *= s;
                                dy *= s;
                                this->positions[2 * particleIndex] -= dx;
                                this->positions[2 * particleIndex + 1] -= dy;
                                this->positions[2 * otherParticleID] += dx;
                                this->positions[2 * otherParticleID + 1] += dy;
                            }
                        }
                    }*/
                }
            }
        }
    }

    void makeForceObjectQueriesConstantMem(bool forceObjectActive) {
        if (forceObjectActive) {
            // might have to make this std::max(1, ...); 
            // using cellSpacing instead of spacing just works for the constant memory hashing idk why
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

    void makeForceObjectQueries(bool forceObjectActive) {
        if (forceObjectActive) {
            // might have to make this std::max(1, ...); 
            float numCovered = int(std::ceil(forceObjectRadius / this->spacing));
            for (int i = -numCovered; i < numCovered + 1; ++i) {
                for (int j = -int(std::min(mouseY / this->spacing, numCovered)); j < int(std::min(this->nY - mouseY / this->spacing, numCovered + 1)); ++j) {
                    int cell = this->getCell(mouseX + i * this->spacing, mouseY + j * this->spacing);

                    if (cell < 0 || cell > tableSize2 - 1) continue;

                    int start = this->cellCount2[cell];
                    int end = this->cellCount2[cell + 1];
                    for (int p = start; p < end; ++p) {
                        int otherParticleID = this->particleArray[p];

                        /*particleColors[3 * otherParticleID] = 0;
                        particleColors[3 * otherParticleID + 1] = 255;
                        particleColors[3 * otherParticleID + 2] = 0;*/

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

    void constrainWalls() {
        for (int i = 0; i < numParticles; ++i) {
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
            }
        }
    }

    void transferVelocities(bool toGrid, float flipRatio) {
        float n = this->numY;
        float h = this->cellSpacing;
        float h1 = this->invSpacing;
        float h2 = 0.5 * h;

        if (toGrid) {
            // make a copy of the grid
            //std::copy(std::begin(this->u), std::end(this->u), std::begin(this->prevU));
            //std::copy(std::begin(this->v), std::end(this->v), std::begin(this->prevV));
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
                int xi = this->clamp(std::floor(x * h1), 0, this->numX - 1);
                int yi = this->clamp(std::floor(y * h1), 0, this->numY - 1);

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

        // now transfer velocities
        for (int component = 0; component < 2; ++component) {
            // u and v grids are staggered, so make sure that you subtract half cell spacing from particle y positions when transferring this->u to particles and vice versa for this->v
            float dx = component == 0 ? 0.f : h2;
            float dy = component == 0 ? h2 : 0.f;

            // on the first pass (component = 0) deal with u, on the second pass (component = 1) deal with v grid
            std::vector<float>* f = component == 0 ? &this->u : &this->v;
           
            // initialize a before grid for the FLIP method
            std::vector<float>* prevF = component == 0 ? &this->prevU : &this->prevV;
           
            std::vector<float>* d = component == 0 ? &this->du : &this->dv;

            for (int i = 0; i < this->numParticles; ++i) {
                float x = this->positions[2 * i];
                float y = this->positions[2 * i + 1];

                x = this->clamp(x, h, (this->numX - 1) * h);
                y = this->clamp(y, h, (this->numY - 1) * h);

                // x0 is the grid position to the left of the particle, x1 is the position to the right of the particle. Both can only go up to the second to last cell to the right in the grid because we dont want to be changing wall velocities
                int x0 = std::max(1, std::min((int)(std::floor((x - dx) * h1)), this->numX - 1));

                // basically x - xCell to get the weight of that cell in relation to the particle
                // in this case, x is moved over to x - dx, and xCell is just grid position of x multiplied by grid spacing
                float tx = ((x - dx) - x0 * h) * h1;

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
                int y0 = std::max(0, std::min((int)(std::floor((y - dy) * h1)), this->numY - 2));
                float ty = ((y - dy) - y0 * h) * h1;
                int y1 = std::min(y0 + 1, this->numY - 1);

                float sx = 1.f - tx;
                float sy = 1.f - ty;

                // weights for each corner in u/v field
                float d0 = sx * sy;
                float d1 = tx * sy;
                float d2 = tx * ty;
                float d3 = sx * ty;

                // top left
                int nr0 = x0 * n + y0;
                // top right
                int nr1 = x1 * n + y0;
                // bottom right
                int nr2 = x1 * n + y1;
                //bottom left
                int nr3 = x0 * n + y1;
               
                if (toGrid) {
                    float pv = this->velocities[2 * i + component];

                    (*f)[nr0] += pv * d0;  
                    (*f)[nr1] += pv * d1;
                    (*f)[nr2] += pv * d2;
                    (*f)[nr3] += pv * d3;

                    (*d)[nr0] += d0;
                    (*d)[nr1] += d1;
                    (*d)[nr2] += d2;
                    (*d)[nr3] += d3;
                }

                else {
                    int offset = component == 0 ? n : 1;
                    // these will be used to make sure that air cells are not considered when transferring velocities back to particles
                    // nr0 - offset is the same this as [(i-1) * n]
                    // if u is being considered, then we only have to check left and right cells ([nr0] and [nr0 - n])
                    // if v is being considered, then we only have to check above and below cells ([nr0] and [nr0 - 1])
                    float valid0 = this->cellType[nr0] != AIR_CELL || this->cellType[nr0 - offset] != AIR_CELL ? 1.0 : 0.0;
                    float valid1 = this->cellType[nr1] != AIR_CELL || this->cellType[nr1 - offset] != AIR_CELL ? 1.0 : 0.0;
                    float valid2 = this->cellType[nr2] != AIR_CELL || this->cellType[nr2 - offset] != AIR_CELL ? 1.0 : 0.0;
                    float valid3 = this->cellType[nr3] != AIR_CELL || this->cellType[nr3 - offset] != AIR_CELL ? 1.0 : 0.0;

                    float v = this->velocities[2 * i + component];
                    float d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                    if (d > 0.0) {
                        float picV = (valid0 * d0 * (*f)[nr0] + valid1 * d1 * (*f)[nr1] + valid2 * d2 * (*f)[nr2] + valid3 * d3 * (*f)[nr3]) / d;

                        float corr = (valid0 * d0 * ((*f)[nr0] - (*prevF)[nr0]) + valid1 * d1 * ((*f)[nr1] - (*prevF)[nr1]) + valid2 * d2 * ((*f)[nr2] - (*prevF)[nr2]) + valid3 * d3 * ((*f)[nr3] - (*prevF)[nr3])) / d;

                        float flipV = v + corr;

                        this->velocities[2 * i + component] = (1.f - flipRatio) * picV + flipRatio * flipV;
                    }
                }
            }

            if (toGrid) {
                for (int i = 0; i < (*f).size(); ++i) {
                    if ((*d)[i] > 0.f) {
                        (*f)[i] /= (*d)[i];
                    }
                }
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


    void solveIncompressibility(const int numIters, const float overRelaxation) {
        std::fill(begin(this->p), end(this->p), 0);

        std::copy(std::begin(this->u), std::end(this->u), std::begin(this->prevU));
        std::copy(std::begin(this->v), std::end(this->v), std::begin(this->prevV));

        int n = this->numY;

        for (int iter = 0; iter < numIters; ++iter) {
            for (int i = 1; i < this->numX - 1; ++i) {
                for (int j = 1; j < this->numY - 1; ++j) {
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
                        float k = 2.3f; // 2.5
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

    // would need a gigantic grid for incompressibility multithreading to be viable 
    /*void solveIncompressibility(const int numIters, const float overRelaxation) {
        std::fill(begin(this->p), end(this->p), 0);

        std::copy(std::begin(this->u), std::end(this->u), std::begin(this->prevU));
        std::copy(std::begin(this->v), std::end(this->v), std::begin(this->prevV));


        for (int iter = 0; iter < numIters; ++iter) {
            this->solvePass(overRelaxation, true);

            this->solvePass(overRelaxation, false);
        }
    }*/

    void includeRigidObject(const bool mouseDown, const bool justPressed, const float dt) {
        if (mouseDown) {
            int n = numY;
            float vx = (objectX - objectPrevX) * 75;
            float vy = (objectY - objectPrevY) * 75;
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

    void updateVertexArray(int index) {
        int i = 4 * index;
        const float px = positions[2 * index];
        const float py = positions[2 * index + 1];

        va[i].position = {px - radius, py - radius};
        va[i + 1].position = {px + radius, py - radius};
        va[i + 2].position = {px + radius, py + radius};
        va[i + 3].position = {px - radius, py + radius};

        sf::Color color;

        int vel = (int)(velocities[2 * index] * velocities[2 * index] + velocities[2 * index + 1] * velocities[2 * index + 1]) / 3000; 
        if (vel > gradient.size()) {
            color = sf::Color(gradient[gradient.size() - 1][0], gradient[gradient.size() - 1][1], gradient[gradient.size() - 1][2], 255);
        }
        else {
            color = sf::Color(gradient[vel][0], gradient[vel][1], gradient[vel][2], 255);
        }

        va[i].color = color;
        va[i + 1].color = color;
        va[i + 2].color = color;
        va[i + 3].color = color;
    }

    void drawParticlesVertex(sf::RenderWindow& window) {
        window.draw(va, states);
    }

    void simulate(float sdt, sf::RenderWindow& window, bool forceObjectActive, bool leftMouseDown, bool justPressed, int numPressureIters, float overRelaxation, float flipRatio, bool rightMouseDown) {

        this->integrate(sdt, window);

        this->initializeSH();
        //this->initializeSHConstantMem();

        for (int i = 0; i < numParticles; ++i) {
            this->makeParticleQueries(i);
        }
        //this->makeParticleQueriesConstantMem(0, numParticles);

        this->makeForceObjectQueries(forceObjectActive);
        //this->makeForceObjectQueriesConstantMem(forceObjectActive);

        this->constrainWalls();

        this->transferVelocities(true, 0.f);

        this->updateParticleDensity();
        if (!forceObjectActive) {
            this->includeRigidObject(leftMouseDown, justPressed, sdt);
        }
        this->solveIncompressibility(numPressureIters, overRelaxation);
        
        this->transferVelocities(false, flipRatio);

        if (forceObjectActive) {
            if (leftMouseDown) {
                this->applyForceObjectForces(250, sdt); // pulling, 250
            }
            else if (rightMouseDown) { // dont worry about the forceObject being stuck in push/pull mode. the boolean in the seperate particles method will take care of that 
                this->applyForceObjectForces(-1000, sdt); // pushing, -1000
            }
        }

        for (int i = 0; i < numParticles; ++i) {
            this->updateVertexArray(i);
        }

        this->drawParticlesVertex(window);

        if (forceObjectActive) {
            this->drawForceObject(window);
        }
        else {
            this->drawRigidObject(window);
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

    void addToForceObjectRadius(float add) {
        if (forceObjectRadius + add > 0) {
            this->forceObjectRadius += add;
            this->forceObjectDrawer.setOrigin(forceObjectRadius, forceObjectRadius);
            this->forceObjectDrawer.setRadius(forceObjectRadius);
            this->checkForceObjectSeperationDist = (this->radius + forceObjectRadius) * (this->radius + forceObjectRadius);
        }
    }
   
};


int main()
{
    int WIDTH = 2000; //2000, 800
    int HEIGHT = 1000; // 800, 1300 
    int numParticles = 10000; // 5000
    float restitution = 0.5f;
    float gravity = 2500.f; //2500

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "FLIP Simulation");

    sf::Font font;
    font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");

    sf::Text text;
    text.setFont(font);
    text.setPosition(10, 10);
    text.setFillColor(sf::Color::White);

    sf::Clock deltaClock;

    window.setFramerateLimit(120);
    window.setMouseCursorVisible(false);

    int frame = 0;
    int fps = 0;

    int subStep = 1;

    float interactionRadius = 100.f;

    float interactionStrength;

    bool leftMouseDown = false;
    bool rightMouseDown = false;

    Fluid fluid = Fluid(1000, WIDTH, HEIGHT, 1.f * HEIGHT / 70, numParticles, gravity);  //50

    float overRelaxation = 1.9f;
    float flipRatio = 0.95f;

    int numPressureIters = 20;

    bool justPressed = false;


    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << flipRatio; 

    bool forceObjectActive = true; 

    float totalDT = 0;
    float numDT = 0;

    while (window.isOpen())
    {
        sf::Time deltaTime = deltaClock.restart();
        float dt = deltaTime.asSeconds();
        float sdt = dt / subStep;

        /*totalDT += dt;
        numDT++;*/

        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::B) {
                    if (flipRatio < 0.99) {
                        flipRatio += 0.01;
                        oss.str("");  
                        oss.clear();
                        oss << std::fixed << std::setprecision(2) << flipRatio; 
                    }
                }
                else if (event.key.code == sf::Keyboard::S) {
                    if (flipRatio > 0.01) {
                        flipRatio -= 0.01;
                        oss.str("");  
                        oss.clear();
                        oss << std::fixed << std::setprecision(2) << flipRatio;
                    }
                }
                else if (event.key.code == sf::Keyboard::Num1) {
                    forceObjectActive = false;
                }
                else if (event.key.code == sf::Keyboard::Num2) {
                    forceObjectActive = true;
                }
                else if (event.key.code == sf::Keyboard::T) {
                    fluid.addToForceObjectRadius(3);
                }
                else if (event.key.code == sf::Keyboard::R) {
                    fluid.addToForceObjectRadius(-3);
                }
                else if (event.key.code == sf::Keyboard::Q) {
                    //std::cout << totalDT / numDT;
                    window.close();
                }
            }
            else if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    leftMouseDown = true;
                    justPressed = true;
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    rightMouseDown = true;
                }
            }
            else if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    leftMouseDown = false;
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    rightMouseDown = false;
                }
            }
        }

        window.clear();

        fluid.simulate(sdt, window, forceObjectActive, leftMouseDown, justPressed, numPressureIters, overRelaxation, flipRatio, rightMouseDown);

        frame++;
        if (frame == 30) {
            fps = (int)(1.f / dt);
            frame = 0;
        }

        text.setPosition(WIDTH - 70, 10);
        text.setString(std::to_string(fps)); 
        window.draw(text);


        text.setPosition(WIDTH - 175, 10);
        text.setString(oss.str());
        window.draw(text);


        window.display();
 
        justPressed = false;
       
    }

    return 0;
}
