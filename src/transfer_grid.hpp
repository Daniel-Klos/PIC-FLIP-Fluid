#pragma once
#include <vector>
#include <cstdint>

#include "fluid_state.hpp"
#include "general_math.hpp"
#include "collision_grid.hpp"

class TransferGrid {

    float invSpacing;
    float halfSpacing;
    int n;

    int SOLID_CELL;
    int FLUID_CELL;
    int AIR_CELL;

    FluidState &fluid_attributes;

    CollisionGrid cellOccupantsGrid;

    std::vector<std::vector<float>> threadDensities;

public:
    std::vector<int32_t> nr0;
    std::vector<int32_t> nr1;
    std::vector<int32_t> nr2;
    std::vector<int32_t> nr3;

    std::vector<float> d0;
    std::vector<float> d1;
    std::vector<float> d2;
    std::vector<float> d3;

    TransferGrid(FluidState &fas): fluid_attributes(fas) {
        invSpacing =  1.f / fluid_attributes.cellSpacing;
        halfSpacing = 0.5f * fluid_attributes.cellSpacing;
        n = fluid_attributes.numY;

        nr0.resize(2 * fluid_attributes.num_particles);
        nr1.resize(2 * fluid_attributes.num_particles);
        nr2.resize(2 * fluid_attributes.num_particles);
        nr3.resize(2 * fluid_attributes.num_particles);
    
        d0.resize(2 * fluid_attributes.num_particles);
        d1.resize(2 * fluid_attributes.num_particles);
        d2.resize(2 * fluid_attributes.num_particles);
        d3.resize(2 * fluid_attributes.num_particles);

        threadDensities.resize(fluid_attributes.numThreads);
        for (int i = 0; i < fluid_attributes.numThreads; ++i) {
            threadDensities[i].resize(fluid_attributes.gridSize, 0.f);
        }

        SOLID_CELL = fluid_attributes.SOLID_CELL;
        FLUID_CELL = fluid_attributes.FLUID_CELL;
        AIR_CELL = fluid_attributes.AIR_CELL;

        cellOccupantsGrid = CollisionGrid(fluid_attributes.numX, fluid_attributes.numY);
    }

    void TransferToGrid() {
        cacheTransferNodesMulti();
        setUpTransferGrids();
        transferParticleVelocitiesToGridMulti();
    }

    void TransferToParticles() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            transferToParticle(start, end);
        });
    }

    void updateCellDensitiesMulti() {
        std::fill(begin(fluid_attributes.cellDensities), end(fluid_attributes.cellDensities), 0.f);
        for (auto& vec : threadDensities) {
            std::fill(vec.begin(), vec.end(), 0.f);
        }

        for (int i = 0; i < fluid_attributes.numThreads; ++i) {
            int start = i * fluid_attributes.particlesPerThread;
            int end = (i == fluid_attributes.numThreads - 1) ? start + fluid_attributes.particlesPerThread + fluid_attributes.numMissedParticles : start + fluid_attributes.particlesPerThread;
            fluid_attributes.thread_pool.addTask([this, start, end, i] {
                updateCellDensities(start, end, i);
            });
        }

        fluid_attributes.thread_pool.waitForCompletion();

        applyLocalDensityUpdatesMulti();

        calculateRestDensity();
    }

private:
    // pgridx = 
    void interpolateTo(float px, float py, float pgridx, float pgridy, float& topL) {

    }

    void cacheWeightsAt(int pIdx, int component, int RKstep) {
        const float dx = (component != 0) * halfSpacing;
        const float dy = (component == 0) * halfSpacing;

        float px = fluid_attributes.positions[2 * pIdx];
        float py = fluid_attributes.positions[2 * pIdx + 1];

        px = clamp(px, fluid_attributes.cellSpacing, (fluid_attributes.numX - 1) * fluid_attributes.cellSpacing);
        py = clamp(py, fluid_attributes.cellSpacing, (fluid_attributes.numY - 1) * fluid_attributes.cellSpacing);

        // x0 is the grid position to the left of the particle, x1 is the position to the right of the particle
        int x0 = std::max(1, std::min(static_cast<int>(std::floor((px - dx) * invSpacing)), fluid_attributes.numX - 2)); // - 1
        // x - xCell to get the lerp weight of that cell in relation to the particle
        float tx = ((px - dx) - x0 * fluid_attributes.cellSpacing) * invSpacing;
        int x1 = std::min(x0 + 1, fluid_attributes.numX - 2); // - 1
        // this fixes a bug that makes water touching the left wall and ceiling explode sometimes -> doesn't allow pos - offset to sample out of bounds
        if (component == 0 && x0 == 1) { // ceiling
            x0 = x1;
        }
        if (component == 1 && x0 == 1) { // wall
            x1 = x0;
        }
        // same thing with y
        int y0 = std::max(0, std::min(static_cast<int>(std::floor((py - dy) * invSpacing)), fluid_attributes.numY - 2)); // - 2
        float ty = ((py - dy) - y0 * fluid_attributes.cellSpacing) * invSpacing;
        int y1 = std::min(y0 + 1, fluid_attributes.numY - 1); // - 1

        float sx = 1.f - tx;
        float sy = 1.f - ty;

        int gIdx = 2 * pIdx + component;
        // weights for each corner in u/v field
        d0[gIdx] = sx * sy; //* fluid_attributes.particle_densities[pIdx];
        d1[gIdx] = tx * sy; //* fluid_attributes.particle_densities[pIdx];
        d2[gIdx] = tx * ty; //* fluid_attributes.particle_densities[pIdx];
        d3[gIdx] = sx * ty; //* fluid_attributes.particle_densities[pIdx];

        nr0[gIdx] = x0 * n + y0; // top left
        nr1[gIdx] = x1 * n + y0; // top right
        nr2[gIdx] = x1 * n + y1; // bottom right
        nr3[gIdx] = x0 * n + y1; // bottom left
    }

    // make nr0 a vector of arrays of size whatever RK you need
    void cacheTransferNodes(int32_t start, int32_t end, float halfHeight, int32_t component) {
        for (int32_t i = start; i < end; ++i) {
            cacheWeightsAt(i, component, 0);
        }
    }

    void cacheTransferNodesMulti() {
        for (int32_t i = 0; i < fluid_attributes.numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([&, i](){
                cacheTransferNodes(fluid_attributes.particlesPerThread * i, fluid_attributes.particlesPerThread * i + fluid_attributes.particlesPerThread, halfSpacing, 0);
                cacheTransferNodes(fluid_attributes.particlesPerThread * i, fluid_attributes.particlesPerThread * i + fluid_attributes.particlesPerThread, halfSpacing, 1);
            });
        }

        cacheTransferNodes(fluid_attributes.num_particles - fluid_attributes.numMissedParticles, fluid_attributes.num_particles, halfSpacing, 0);
        cacheTransferNodes(fluid_attributes.num_particles - fluid_attributes.numMissedParticles, fluid_attributes.num_particles, halfSpacing, 1);

        fluid_attributes.thread_pool.waitForCompletion();
    }

    void setUpTransferGrids() {
        const float h2 = 0.5 * fluid_attributes.cellSpacing;

        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));
        std::fill(begin(fluid_attributes.du), end(fluid_attributes.du), 0.f);
        std::fill(begin(fluid_attributes.dv), end(fluid_attributes.dv), 0.f);
        std::fill(begin(fluid_attributes.u), end(fluid_attributes.u), 0.f);
        std::fill(begin(fluid_attributes.v), end(fluid_attributes.v), 0.f);

        // initialize every inside cell to air
        for (int i = 0; i < fluid_attributes.numX; ++i) {
            for (int j = 0; j < fluid_attributes.numY; ++j) {
                if (fluid_attributes.cellType[i * n + j] != SOLID_CELL) {
                    fluid_attributes.cellType[i * n + j] = AIR_CELL;
                }
            }
        }

        // initialize all cells that particles are in to fluid
        fluid_attributes.num_fluid_cells = 0;
        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            float x = fluid_attributes.positions[2 * i];
            float y = fluid_attributes.positions[2 * i + 1];

            int xi = clamp(std::floor(x * invSpacing), 0, fluid_attributes.numX - 1);
            int yi = clamp(std::floor(y * invSpacing), 0, fluid_attributes.numY - 1);

            int cellNr = xi * n + yi;
            if (fluid_attributes.cellType[cellNr] == AIR_CELL) {
                fluid_attributes.cellType[cellNr] = FLUID_CELL;
                fluid_attributes.fluidCellPositions[fluid_attributes.num_fluid_cells] = sf::Vector2i{xi, yi};
                fluid_attributes.num_fluid_cells++;
            }
        }

        // add particles to occupantsGrid for multithreading
        cellOccupantsGrid.clear();

        const float minX = fluid_attributes.cellSpacing;
        const float maxX = fluid_attributes.WIDTH - fluid_attributes.cellSpacing;
        const float minY = fluid_attributes.cellSpacing;
        const float maxY = fluid_attributes.HEIGHT - fluid_attributes.cellSpacing;

        uint32_t i{0};

        for (int32_t index = 0; index < fluid_attributes.num_particles; ++index) {
            float x = fluid_attributes.positions[2 * index];
            float y = fluid_attributes.positions[2 * index + 1];
            if (x > minX && x < maxX && y > minY && y < maxY) {
                
                int32_t cellOccupantsX = x / fluid_attributes.cellSpacing;
                int32_t cellOccupantsY = y / fluid_attributes.cellSpacing;
                cellOccupantsGrid.addAtom(cellOccupantsX, cellOccupantsY, i);

            }
            ++i;
        }
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
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
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
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int idx = i * n + j;
                transferToVGridCells(idx);
            }
        }
    }

    void normalizeGridInterpolations(int start, int end) {
        for (int i = start; i < end; ++i) {
            float prevNode = fluid_attributes.du[i];
            if (prevNode > 0.f) {
                fluid_attributes.u[i] /= prevNode;
            }
            prevNode = fluid_attributes.dv[i];
            if (prevNode > 0.f) {
                fluid_attributes.v[i] /= prevNode;
            }
        }
    }

    void normalizeGridInterpolationsMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.gridSize, [this](int start, int end) {
            normalizeGridInterpolations(start, end);
        });
    }

    void enforceNoSlip(int start, int end) {
        for (int i = start; i < end; ++i) {
            for (int j = 0; j < fluid_attributes.numY; ++j) {
                int idx = i * n + j;
                bool solid = fluid_attributes.cellType[idx] == SOLID_CELL;
                if (solid || i > 0 && fluid_attributes.cellType[idx - n] == SOLID_CELL) {
                    fluid_attributes.u[idx] = 0;
                }
                if (solid || j > 0 && fluid_attributes.cellType[idx - 1] == SOLID_CELL) {
                    fluid_attributes.v[idx] = 0;
                }
            }
        }
    }

    void enforceNoSlipMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.numX, [this](int start, int end) {
            enforceNoSlip(start, end);
        });
    }

    void transferParticleVelocitiesToGridMulti() {

        const int32_t numThreads = fluid_attributes.numThreads;
        const int32_t numColumnsPerThread = (fluid_attributes.numX - 2) / numThreads;
        const int32_t numMissedColumns = fluid_attributes.numX - 2 - numColumnsPerThread * numThreads;

        for (int i = 0; i < numThreads; ++i) {
            int start = i * numColumnsPerThread;
            int end = (i == numThreads - 1) ? (fluid_attributes.numX - 1) : (start + numColumnsPerThread);
            fluid_attributes.thread_pool.addTask([&, start, end]() {
                transferToUGrid(start, end);
            });
        }

        for (int i = 0; i < numThreads; ++i) {
            int start = i * numColumnsPerThread;
            int end = (i == numThreads - 1) ? (fluid_attributes.numX - 1) : (start + numColumnsPerThread);
            fluid_attributes.thread_pool.addTask([&, start, end]() {
                transferToVGrid(start, end);
            });
        }

        fluid_attributes.thread_pool.waitForCompletion();

        normalizeGridInterpolationsMulti();

        enforceNoSlipMulti();

        // store prev grid for FLIP
        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));
    }



    void transferToParticle(int32_t startIndex, int32_t endIndex) {
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
           
            // AIR_CELL check logic
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

    void calculateRestDensity() {
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

    void updateCellDensities(int start, int end, int threadID) {
        for (int i = start; i < end; ++i) {
            float x = fluid_attributes.positions[2 * i];
            float y = fluid_attributes.positions[2 * i + 1];

            x = clamp(x, fluid_attributes.cellSpacing, (fluid_attributes.numX - 1) * fluid_attributes.cellSpacing);
            y = clamp(y, fluid_attributes.cellSpacing, (fluid_attributes.numY - 1) * fluid_attributes.cellSpacing);

            int x0 = std::max(1, std::min((int)(std::floor((x - halfSpacing) * invSpacing)), fluid_attributes.numX - 2));
            float tx = ((x - halfSpacing) - x0 * fluid_attributes.cellSpacing) * invSpacing;
            int x1 = std::min(x0 + 1, fluid_attributes.numX - 1);

            int y0 = std::max(1, std::min((int)(std::floor((y - halfSpacing) * invSpacing)), fluid_attributes.numY - 2));
            float ty = ((y - halfSpacing) - y0 * fluid_attributes.cellSpacing) * invSpacing;
            int y1 = std::min(y0 + 1, fluid_attributes.numY - 2);
           
            float sx = 1.f - tx;
            float sy = 1.f - ty;

            if (x0 < fluid_attributes.numX && y0 < fluid_attributes.numY) {
                threadDensities[threadID][x0 * n + y0] += sx * sy;
            }
            if (x1 < fluid_attributes.numX && y0 < fluid_attributes.numY) {
                threadDensities[threadID][x1 * n + y0] += tx * sy;
            }
            if (x1 < fluid_attributes.numX && y1 < fluid_attributes.numY) {
                threadDensities[threadID][x1 * n + y1] += tx * ty;
            }
            if (x0 < fluid_attributes.numX && y1 < fluid_attributes.numY) {
                threadDensities[threadID][x0 * n + y1] += sx * ty;
            }
        }
    }

    void applyLocalDensityUpdates(int start, int end) {
        for (int i = start; i < end; ++i) {
            for (int j = 0; j < fluid_attributes.numThreads; ++j) {
                fluid_attributes.cellDensities[i] += threadDensities[j][i];
            }
        }
    }

    void applyLocalDensityUpdatesMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.gridSize, [this](int start, int end) {
            applyLocalDensityUpdates(start, end);
        });
    }

    /*void updateParticleDensities(int start, int end) {
        for (int i = start; i < end; ++i) {
            float px = clamp(fluid_attributes.positions[2 * i], fluid_attributes.cellSpacing, (fluid_attributes.numX - 1) * fluid_attributes.cellSpacing);
            float py = clamp(fluid_attributes.positions[2 * i + 1], fluid_attributes.cellSpacing, (fluid_attributes.numY - 1) * fluid_attributes.cellSpacing);

            int x0 = std::max(1, std::min((int)std::floor(px * invSpacing), fluid_attributes.numX - 2));
            int y0 = std::max(1, std::min((int)std::floor(py * invSpacing), fluid_attributes.numY - 2));
            int x1 = x0 + 1;
            int y1 = y0 + 1;

            int topL = x0 * n + y0;
            int topR = x1 * n + y0;
            int botL = x0 * n + y1;
            int botR = x1 * n + y1;

            float tx = (px - x0 * fluid_attributes.cellSpacing) * invSpacing;
            float ty = (py - y0 * fluid_attributes.cellSpacing) * invSpacing;
            float sx = 1.f - tx;
            float sy = 1.f - ty;

            fluid_attributes.densities[i] = 0.f;
            fluid_attributes.densities[i] += sx * sy * fluid_attributes.cellDensities[topL];
            fluid_attributes.densities[i] += tx * sy * fluid_attributes.cellDensities[topR];
            fluid_attributes.densities[i] += sx * ty * fluid_attributes.cellDensities[botL];
            fluid_attributes.densities[i] += tx * ty * fluid_attributes.cellDensities[botR];
        }
    }*/
};