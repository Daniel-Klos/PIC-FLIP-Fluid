#pragma once
#include <vector>
#include <cmath>

#include "thread_pool.hpp"

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
    std::vector<float> cellDensities;
    float cellSpacing;
    float radius;
    int numX;
    int numY;
    float WIDTH;
    float HEIGHT;
    int gridSize;

    float particleRestDensity;

    float gravityX;
    float gravityY;

    float vorticityStrength;
    float flipRatio; // move this to TransferGrid struct

    bool stop = false;
    bool step = false; 

    bool fireActive = false;
    std::vector<float> temperatures;
    float groundConductivity = 30000.f;  // how quickly the ground transfers heat to the particles
    float interConductivity = 10000.f;    // how quickly particles transfer heat between themselves
    float fireStrength = 175.f;         // how quickly particles accelerate upwards due to heat
    float tempDiffusion = 0.1f;        // how quickly particles lose heat

    int FLUID_CELL = 0;
    int AIR_CELL = 1;
    int SOLID_CELL = 2;

    tp::ThreadPool& thread_pool;
    int numThreads;
    int particlesPerThread;
    int numMissedParticles;

    FluidState(int num_particles_, float WIDTH_, float HEIGHT_, int numX_, float vorticityStrength_, float flipRatio_, float gravityX_, float gravityY_, tp::ThreadPool &tp): num_particles(num_particles_), WIDTH(WIDTH_), HEIGHT(HEIGHT_), numX(numX_), vorticityStrength(vorticityStrength_), flipRatio(flipRatio_), gravityX(gravityX_), gravityY(gravityY_), thread_pool(tp) {

        cellSpacing = WIDTH / numX;
        radius = 0.3 * cellSpacing;

        numY = std::floor(HEIGHT / cellSpacing);
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

        numThreads = thread_pool.m_thread_count;
        particlesPerThread = num_particles / numThreads;
        numMissedParticles = num_particles - particlesPerThread * numThreads;

        // initializing particle positions
        float separation = 2.1; 

        int wideNum = std::floor((WIDTH - 2 * cellSpacing - 2) / (radius * 2.1));
        int highNum = num_particles / wideNum;

        float starting_px = radius + cellSpacing + 2;
        float starting_py = (HEIGHT - (radius * separation * highNum)) / 2 + radius;

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
};
