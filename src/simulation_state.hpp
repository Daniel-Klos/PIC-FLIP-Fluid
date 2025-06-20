#pragma once
#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>

#include "thread_pool.hpp"
#include "frame_context.hpp"

struct FluidState {
    int num_particles;
    std::vector<float> positions;
    std::vector<float> velocities;
    std::vector<int> cellType;
    std::vector<float> u;
    std::vector<float> v;
    std::vector<float> du;
    std::vector<float> dv;
    std::vector<float> prevU;
    std::vector<float> prevV;
    std::vector<float> cellDensities;       // density at each cell center
    std::vector<int> fluid_cells;
    std::vector<sf::Vector2i> obstaclePositions;
    //std::vector<std::array<float, 4>> Cp;
    //std::vector<float> densities;           // density at each particle
    //std::vector<float> particle_densities;  // density of each particle
    float cellSpacing;
    float halfSpacing;
    float invSpacing;
    float radius;
    int numX;
    int numY;
    int n; // same as numY
    int gridSize;
    int num_fluid_cells;

    float particleRestDensity;

    float gravityX;
    float gravityY;

    float vorticityStrength;
    float flipRatio; // move this to TransferGrid struct

    bool stop = false;
    bool step = false; 

    bool fireActive = false;
    std::vector<float> temperatures;
    float groundConductivity = 20000.f;  // how quickly the ground transfers heat to the particles
    float interConductivity = 10000.f;    // how quickly particles transfer heat between themselves
    float fireStrength = 75.f;         // how quickly particles accelerate upwards due to heat
    float tempDiffusion = 50.f;        // how quickly particles lose heat

    int FLUID_CELL = 0;
    int AIR_CELL = 1;
    int SOLID_CELL = 2;

    tp::ThreadPool& thread_pool;
    int numThreads;
    int particlesPerThread;
    int numMissedParticles;

    FrameContext &frame_context;

    FluidState(int num_particles_, int numX_, float vorticityStrength_, float flipRatio_, float gravityX_, float gravityY_, FrameContext &fc, tp::ThreadPool &tp): num_particles(num_particles_), numX(numX_), vorticityStrength(vorticityStrength_), flipRatio(flipRatio_), gravityX(gravityX_), gravityY(gravityY_), frame_context(fc), thread_pool(tp) {


        cellSpacing = frame_context.WIDTH / numX;
        halfSpacing = cellSpacing * 0.5;
        invSpacing = 1.f / cellSpacing;
        radius = 0.3 * cellSpacing;

        numY = std::floor(frame_context.HEIGHT / cellSpacing);
        n = numY;
        gridSize = numX * numY;

        positions.resize(2 * num_particles);
        velocities.resize(2 * num_particles);
        temperatures.resize(num_particles);
        cellType.resize(gridSize);
        u.resize(gridSize);
        v.resize(gridSize);
        du.resize(gridSize);
        dv.resize(gridSize);
        prevU.resize(gridSize);
        prevV.resize(gridSize);
        cellDensities.resize(gridSize);
        fluid_cells.resize(gridSize);
        /*densities.resize(num_particles);
        particle_densities.resize(num_particles);*/

        numThreads = thread_pool.m_thread_count;
        particlesPerThread = num_particles / numThreads;
        numMissedParticles = num_particles - particlesPerThread * numThreads;

        // initialize all attributes you need to initialize for particles before sim starts
        /*for (int i = 0; i < num_particles; ++i) {
            particle_densities[i] = (i < num_particles / 2) ? 0.9f : 0.1f;
        }*/

        // initializing particle positions
        float separation = 2.1; 

        int wideNum = std::floor((frame_context.WIDTH - 2 * cellSpacing - 2) / (radius * 2.1));
        int highNum = num_particles / wideNum;

        float starting_px = radius + cellSpacing + 2;
        float starting_py = (frame_context.HEIGHT - (radius * separation * highNum)) / 2 + radius;

        float px = starting_px;
        float py = starting_py;

        int addTo = num_particles - (wideNum * highNum);

        bool offset = true;
        for (int i = 0; i < wideNum * highNum + addTo; ++i) {
            this->positions[i * 2] = px;
            this->positions[i * 2 + 1] = py;

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
    }

    float curl(int i, int j) {
        int idx = i * n + j;
        const float denom = 1.f / (2.f * cellSpacing);
        const int32_t leftType = cellType[idx - n] == 0;
        const int32_t rightType = cellType[idx + n] == 0;
        const int32_t topType = cellType[idx - 1] == 0;
        const int32_t bottomType = cellType[idx + 1] == 0;
        if (!leftType || !rightType || !topType || !bottomType) {
            return 0.f;
        }
        return ((v[(i + 1) * n + j] * bottomType - v[(i - 1) * n + j] * topType) - (u[i * n + j + 1] * rightType - u[i * n + j - 1] * leftType)) * denom;
    }

    float calcVorticity(int i, int j) {
        float curl = this->curl(i, j);
        return curl * curl; // std::abs for just curl
    }

    sf::Vector2f gridCellToPos(int idx) {
        int i = idx % numY;
        int j = idx / numY;
        float x = halfSpacing + (j * cellSpacing);
        float y = halfSpacing + (i * cellSpacing);
        return sf::Vector2f{x, y};
    }

    int getNumParticles() {
        return this->num_particles;
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

    void addToVorticityStrength(float add) {
        this->vorticityStrength += add;
    }

    float getVorticityStrength() {
        return this->vorticityStrength;
    }

    float getFlipRatio() {
        return this->flipRatio;
    }

    void addToFlipRatio(const float add) {
        this->flipRatio += add;
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

    void setFireActive(bool active) {
        this->fireActive = active;
    }

    float getFireActive() {
        return this->fireActive;
    }

    int getNumX() {
        return this->numX;
    }

    int getNumY() {
        return this->numY;
    }
};