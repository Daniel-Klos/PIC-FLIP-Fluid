#pragma once
#include <vector>
#include <cstdint>

#include "simulation_state.hpp"
#include "utils.hpp"
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

    std::vector<std::vector<float>> threadGridBuffers;

public:
    std::vector<int32_t> topLeftCells;
    std::vector<int32_t> topRightCells;
    std::vector<int32_t> bottomRightCells;
    std::vector<int32_t> bottomLeftCells;

    std::vector<float> topLeftWeights;
    std::vector<float> topRightWeights;
    std::vector<float> bottomLeftWeights;
    std::vector<float> bottomRightWeights;

    TransferGrid(FluidState &fas): fluid_attributes(fas) {
        invSpacing =  1.f / fluid_attributes.cellSpacing;
        halfSpacing = 0.5f * fluid_attributes.cellSpacing;
        n = fluid_attributes.numY;

        topLeftCells.resize(2 * fluid_attributes.num_particles);
        topRightCells.resize(2 * fluid_attributes.num_particles);
        bottomRightCells.resize(2 * fluid_attributes.num_particles);
        bottomLeftCells.resize(2 * fluid_attributes.num_particles);
    
        topLeftWeights.resize(2 * fluid_attributes.num_particles);
        topRightWeights.resize(2 * fluid_attributes.num_particles);
        bottomLeftWeights.resize(2 * fluid_attributes.num_particles);
        bottomRightWeights.resize(2 * fluid_attributes.num_particles);

        threadGridBuffers.resize(fluid_attributes.numThreads);
        for (int i = 0; i < fluid_attributes.numThreads; ++i) {
            threadGridBuffers[i].resize(fluid_attributes.gridSize, 0.f);
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
            transferUGridToParticles(start, end);
            transferVGridToParticles(start, end);
        });
    }

    void updateCellDensitiesMulti() {
        std::fill(begin(fluid_attributes.cellDensities), end(fluid_attributes.cellDensities), 0.f);
        for (auto& vec : threadGridBuffers) {
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

    float Weight(float distX, float distY) {
        return distX * distY;
    }

    float WeightDerivative(bool positive) {
        return -1.f + 2.f * positive;
    }

    // add RKstep as a parameter
    void cacheWeightsAt(int pIdx, int component) {
        const float dx = (component == 1) * halfSpacing;
        const float dy = (component == 0) * halfSpacing;

        float px = fluid_attributes.positions[2 * pIdx];
        float py = fluid_attributes.positions[2 * pIdx + 1];

        px = clamp(px, fluid_attributes.cellSpacing, (fluid_attributes.numX - 1) * fluid_attributes.cellSpacing);
        py = clamp(py, fluid_attributes.cellSpacing, (fluid_attributes.numY - 1) * fluid_attributes.cellSpacing);

        // set first argument to 0 if you want fluid touching left wall to go down
        int x0 = std::max(1, std::min(static_cast<int>(std::floor((px - dx) * invSpacing)), fluid_attributes.numX - 2));
        int x1 = std::min(x0 + 1, fluid_attributes.numX - 2);

        // this stops the fluid touching the left wall from going down
        if (component == 0 && x0 == 1) { // ceiling
            x0 = x1;
        }
        if (component == 1 && x0 == 1) { // left wall
            x1 = x0;
        }
        
        int y0 = std::max(0, std::min(static_cast<int>(std::floor((py - dy) * invSpacing)), fluid_attributes.numY - 2));
        int y1 = std::min(y0 + 1, fluid_attributes.numY - 1);

        float dxLeft = ((px - dx) - x0 * fluid_attributes.cellSpacing) * invSpacing;
        float dyTop = ((py - dy) - y0 * fluid_attributes.cellSpacing) * invSpacing;
        float dxRight = 1.f - dxLeft;
        float dyBottom = 1.f - dyTop;

        fluid_attributes.dxLefts[pIdx].vu[component] = -dxLeft;
        fluid_attributes.dyTops[pIdx].vu[component] = -dyTop;
        fluid_attributes.dxRights[pIdx].vu[component] = dxRight;
        fluid_attributes.dyBottoms[pIdx].vu[component] = dyBottom;

        int gIdx = 2 * pIdx + component;
    
        topLeftWeights[gIdx] = Weight(dxRight, dyBottom); //* fluid_attributes.masses[pIdx];
        topRightWeights[gIdx] = Weight(dxLeft, dyBottom); //* fluid_attributes.masses[pIdx];
        bottomLeftWeights[gIdx] = Weight(dxRight, dyTop); //* fluid_attributes.masses[pIdx];
        bottomRightWeights[gIdx] = Weight(dxLeft, dyTop); //* fluid_attributes.masses[pIdx];

        topLeftCells[gIdx] = x0 * n + y0;
        topRightCells[gIdx] = x1 * n + y0;
        bottomLeftCells[gIdx] = x0 * n + y1;
        bottomRightCells[gIdx] = x1 * n + y1;
    }

    // make topLeftCell a vector of arrays of size whatever RK you need
    void cacheTransferNodes(int32_t start, int32_t end, int32_t component) {
        for (int32_t i = start; i < end; ++i) {
            cacheWeightsAt(i, component);
        }
    }

    void cacheTransferNodesMulti() {
        for (int32_t i = 0; i < fluid_attributes.numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([&, i](){
                cacheTransferNodes(fluid_attributes.particlesPerThread * i, fluid_attributes.particlesPerThread * i + fluid_attributes.particlesPerThread, 0);
                cacheTransferNodes(fluid_attributes.particlesPerThread * i, fluid_attributes.particlesPerThread * i + fluid_attributes.particlesPerThread, 1);
            });
        }

        cacheTransferNodes(fluid_attributes.num_particles - fluid_attributes.numMissedParticles, fluid_attributes.num_particles, 0);
        cacheTransferNodes(fluid_attributes.num_particles - fluid_attributes.numMissedParticles, fluid_attributes.num_particles, 1);

        fluid_attributes.thread_pool.waitForCompletion();
    }

    void setUpTransferGrids() {
        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));
        std::fill(begin(fluid_attributes.uGridWeights), end(fluid_attributes.uGridWeights), 0.f);
        std::fill(begin(fluid_attributes.vGridWeights), end(fluid_attributes.vGridWeights), 0.f);
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

        // initialize all cells that particles are in to fluid and store fluid cells
        fluid_attributes.num_fluid_cells = 0;
        for (int i = 0; i < fluid_attributes.num_particles; ++i) {
            float x = fluid_attributes.positions[2 * i];
            float y = fluid_attributes.positions[2 * i + 1];

            int xi = clamp(std::floor(x * invSpacing), 0, fluid_attributes.numX - 1);
            int yi = clamp(std::floor(y * invSpacing), 0, fluid_attributes.numY - 1);

            int cellNr = xi * n + yi;
            if (fluid_attributes.cellType[cellNr] == AIR_CELL) {
                fluid_attributes.cellType[cellNr] = FLUID_CELL;
                fluid_attributes.fluid_cells[fluid_attributes.num_fluid_cells] = cellNr;
                fluid_attributes.num_fluid_cells++;
            }
        }

        // add particles to occupantsGrid for multithreading
        cellOccupantsGrid.clear();

        const float minX = fluid_attributes.cellSpacing;
        const float maxX = fluid_attributes.frame_context.WIDTH - fluid_attributes.cellSpacing;
        const float minY = fluid_attributes.cellSpacing;
        const float maxY = fluid_attributes.frame_context.HEIGHT - fluid_attributes.cellSpacing;

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
            const uint32_t pIdx = cell.objects[id];

            const int32_t ui = 2 * pIdx;

            int32_t topLeftCell_u = topLeftCells[ui];
            int32_t topRightCell_u = topRightCells[ui];
            int32_t bottomRightCell_u = bottomRightCells[ui];
            int32_t bottomLeftCell_u = bottomLeftCells[ui];

            float topLeftWeight_u = topLeftWeights[ui];
            float topRightWeight_u = topRightWeights[ui];
            float bottomLeftWeight_u = bottomLeftWeights[ui];
            float bottomRightWeight_u = bottomRightWeights[ui];

            float pvx = fluid_attributes.velocities[ui];

            float leftDist = fluid_attributes.dxLefts[pIdx].vu[1];
            float rightDist = fluid_attributes.dxRights[pIdx].vu[1];
            float bottomDist = fluid_attributes.dyBottoms[pIdx].vu[1];
            float topDist = fluid_attributes.dyTops[pIdx].vu[1];

            int matIdx = 4 * pIdx;

            float C11 = fluid_attributes.affineMats[matIdx];
            float C12 = fluid_attributes.affineMats[matIdx + 1];

            float affinepvxBottomLeft = pvx + 
                                    C11 * leftDist + 
                                    C12 * bottomDist;

            float affinepvxBottomRight = pvx + 
                                    C11 * rightDist + 
                                    C12 * bottomDist;

            float affinepvxTopLeft = pvx + 
                                    C11 * leftDist + 
                                    C12 * topDist;

            float affinepvTopRight = pvx + 
                                    C11 * rightDist + 
                                    C12 * topDist;
           
            /*fluid_attributes.u[topLeftCell_u] += affinepvxTopLeft * topLeftWeight_u;  
            fluid_attributes.u[topRightCell_u] += affinepvTopRight * topRightWeight_u;
            fluid_attributes.u[bottomLeftCell_u] += affinepvxBottomLeft * bottomLeftWeight_u;
            fluid_attributes.u[bottomRightCell_u] += affinepvxBottomRight * bottomRightWeight_u;*/

            fluid_attributes.u[topLeftCell_u] += pvx * topLeftWeight_u;  
            fluid_attributes.u[topRightCell_u] += pvx * topRightWeight_u;
            fluid_attributes.u[bottomLeftCell_u] += pvx * bottomLeftWeight_u;
            fluid_attributes.u[bottomRightCell_u] += pvx * bottomRightWeight_u;

            fluid_attributes.uGridWeights[topLeftCell_u] += topLeftWeight_u;
            fluid_attributes.uGridWeights[topRightCell_u] += topRightWeight_u;
            fluid_attributes.uGridWeights[bottomLeftCell_u] += bottomLeftWeight_u;
            fluid_attributes.uGridWeights[bottomRightCell_u] += bottomRightWeight_u;
        }
    }

    void transferToUGrid(int start, int end) {
        for (int i = start; i < end; ++i) {
            transferToUGridCells(fluid_attributes.fluid_cells[i]);
        }
    }

    void transferToVGridCells(int idx) {
        const auto cell = cellOccupantsGrid.data[idx];
    
        for (uint32_t id{0}; id < cell.objects_count; ++id) {
            const uint32_t pIdx = cell.objects[id];

            const int32_t vi = 2 * pIdx + 1;

            const int32_t topLeftCell_v = topLeftCells[vi];
            const int32_t topRightCell_v = topRightCells[vi];
            const int32_t bottomRightCell_v = bottomRightCells[vi];
            const int32_t bottomLeftCell_v = bottomLeftCells[vi];

            const float topLeftWeight_v = topLeftWeights[vi];
            const float topRightWeight_v = topRightWeights[vi];
            const float bottomLeftWeight_v = bottomLeftWeights[vi];
            const float bottomRightWeight_v = bottomRightWeights[vi];

            const float pvy = fluid_attributes.velocities[vi];

            float leftDist = fluid_attributes.dxLefts[pIdx].vu[0];
            float rightDist = fluid_attributes.dxRights[pIdx].vu[0];
            float bottomDist = fluid_attributes.dyBottoms[pIdx].vu[0];
            float topDist = fluid_attributes.dyTops[pIdx].vu[0];

            int matIdx = 4 * pIdx;

            float C21 = fluid_attributes.affineMats[matIdx + 2];
            float C22 = fluid_attributes.affineMats[matIdx + 3];

            float affinepvyBottomLeft = pvy + 
                                    C21 * leftDist + 
                                    C22 * bottomDist;

            float affinepvyBottomRight = pvy + 
                                    C21 * rightDist + 
                                    C22 * bottomDist;

            float affinepvyTopLeft = pvy + 
                                    C21 * leftDist + 
                                    C22 * topDist;

            float affinepvyTopRight = pvy + 
                                    C21 * rightDist + 
                                    C22 * topDist;
    
            /*fluid_attributes.v[topLeftCell_v] += affinepvyTopLeft * topLeftWeight_v;  
            fluid_attributes.v[topRightCell_v] += affinepvyTopRight * topRightWeight_v;
            fluid_attributes.v[bottomLeftCell_v] += affinepvyBottomLeft * bottomLeftWeight_v;
            fluid_attributes.v[bottomRightCell_v] += affinepvyBottomRight * bottomRightWeight_v;*/

            fluid_attributes.v[topLeftCell_v] += pvy * topLeftWeight_v;  
            fluid_attributes.v[topRightCell_v] += pvy * topRightWeight_v;
            fluid_attributes.v[bottomLeftCell_v] += pvy * bottomLeftWeight_v;
            fluid_attributes.v[bottomRightCell_v] += pvy * bottomRightWeight_v;

            fluid_attributes.vGridWeights[topLeftCell_v] += topLeftWeight_v;
            fluid_attributes.vGridWeights[topRightCell_v] += topRightWeight_v;
            fluid_attributes.vGridWeights[bottomLeftCell_v] += bottomLeftWeight_v;
            fluid_attributes.vGridWeights[bottomRightCell_v] += bottomRightWeight_v;
        }
    }

    void transferToVGrid(int start, int end) {
        for (int i = start; i < end; ++i) {
            transferToVGridCells(fluid_attributes.fluid_cells[i]);
        }
    }

    void momentumToVelocity(int start, int end) {
        for (int i = start; i < end; ++i) {
            float prevNode = fluid_attributes.uGridWeights[i];
            if (prevNode > 0.f) {
                fluid_attributes.u[i] /= prevNode;
            }
            prevNode = fluid_attributes.vGridWeights[i];
            if (prevNode > 0.f) {
                fluid_attributes.v[i] /= prevNode;
            }
        }
    }

    void momentumToVelocityMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.gridSize, [this](int start, int end) {
            momentumToVelocity(start, end);
        });
    }

    /*
    void normalizeGridInterpolations(int start, int end) {
        for (int i = start; i < end; ++i) {
            int cellNr = fluid_attributes.fluid_cells[i];
            
            float prevNode = fluid_attributes.du[cellNr];
            if (prevNode > 0.f) {
                fluid_attributes.u[cellNr] /= prevNode;
            }
            prevNode = fluid_attributes.dv[cellNr];
            if (prevNode > 0.f) {
                fluid_attributes.v[cellNr] /= prevNode;
            }
        }
    }

    void normalizeGridInterpolationsMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_fluid_cells, [this](int start, int end) {
            normalizeGridInterpolations(start, end);
        });
    }
    */

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
        const int32_t numCellsPerThread = fluid_attributes.num_fluid_cells / numThreads;
        const int32_t numMissedColumns = fluid_attributes.num_fluid_cells - numCellsPerThread * numThreads;

        for (int i = 0; i < numThreads; ++i) {
            int start = i * numCellsPerThread;
            int end = (i == numThreads - 1) ? (fluid_attributes.num_fluid_cells) : (start + numCellsPerThread);
            fluid_attributes.thread_pool.addTask([&, start, end]() {
                transferToUGrid(start, end);
            });
        }

        for (int i = 0; i < numThreads; ++i) {
            int start = i * numCellsPerThread;
            int end = (i == numThreads - 1) ? (fluid_attributes.num_fluid_cells) : (start + numCellsPerThread);
            fluid_attributes.thread_pool.addTask([&, start, end]() {
                transferToVGrid(start, end);
            });
        }

        fluid_attributes.thread_pool.waitForCompletion();

        momentumToVelocityMulti();

        enforceNoSlipMulti();

        // store prev grid for FLIP
        std::copy(std::begin(fluid_attributes.u), std::end(fluid_attributes.u), std::begin(fluid_attributes.prevU));
        std::copy(std::begin(fluid_attributes.v), std::end(fluid_attributes.v), std::begin(fluid_attributes.prevV));
    }

    void transferUGridToParticles(int32_t startIndex, int32_t endIndex) {
        for (int32_t i = startIndex; i < endIndex; ++i) {
            int ui = 2 * i;
            float pvx = fluid_attributes.velocities[ui];
            int matIdx = 4 * i;

            int topLeftCell = topLeftCells[ui];
            int topRightCell = topRightCells[ui];
            int bottomRightCell = bottomRightCells[ui];
            int bottomLeftCell = bottomLeftCells[ui];

            float topLeftWeight = topLeftWeights[ui];
            float topRightWeight = topRightWeights[ui];
            float bottomLeftWeight = bottomLeftWeights[ui];
            float bottomRightWeight = bottomRightWeights[ui];

            float gridVel_topLeft = fluid_attributes.u[topLeftCell];
            float gridVel_topRight = fluid_attributes.u[topRightCell];
            float gridVel_bottomRight = fluid_attributes.u[bottomRightCell];
            float gridVel_bottomLeft = fluid_attributes.u[bottomLeftCell];
           
            const float validTopLeftWeight = fluid_attributes.cellType[topLeftCell] != AIR_CELL || 
                                             fluid_attributes.cellType[topLeftCell - n] != AIR_CELL;
            const float validTopRightWeight = fluid_attributes.cellType[topRightCell] != AIR_CELL || 
                                              fluid_attributes.cellType[topRightCell - n] != AIR_CELL;
            const float validBottomLeftWeight = fluid_attributes.cellType[bottomLeftCell] != AIR_CELL || 
                                                fluid_attributes.cellType[bottomLeftCell - n] != AIR_CELL;
            const float validBottomRightWeight = fluid_attributes.cellType[bottomRightCell] != AIR_CELL || 
                                                 fluid_attributes.cellType[bottomRightCell - n] != AIR_CELL;

            const float divX = validTopLeftWeight * topLeftWeight + validTopRightWeight * topRightWeight + validBottomLeftWeight * bottomLeftWeight + validBottomRightWeight * bottomRightWeight;

            float picV;
            float corr;
            float flipV;

            if (divX > 0.f) {
                picV = (validTopLeftWeight * topLeftWeight * gridVel_topLeft + 
                        validTopRightWeight * topRightWeight * gridVel_topRight + 
                        validBottomLeftWeight * bottomLeftWeight * gridVel_bottomLeft + 
                        validBottomRightWeight * bottomRightWeight * gridVel_bottomRight
                        ) / divX;

                corr = (validTopLeftWeight * topLeftWeight * (gridVel_topLeft - fluid_attributes.prevU[topLeftCell]) + 
                        validTopRightWeight * topRightWeight * (gridVel_topRight - fluid_attributes.prevU[topRightCell]) + 
                        validBottomLeftWeight * bottomLeftWeight * (gridVel_bottomLeft - fluid_attributes.prevU[bottomLeftCell]) + 
                        validBottomRightWeight * bottomRightWeight * (gridVel_bottomRight - fluid_attributes.prevU[bottomRightCell])
                        ) / divX;

                flipV = pvx + corr;
                fluid_attributes.velocities[ui] = (1.f - fluid_attributes.flipRatio) * picV + fluid_attributes.flipRatio * flipV;
            }

            float weightGradXLeft = fluid_attributes.dxLefts[i].vu[1];
            float weightGradYTop = fluid_attributes.dyTops[i].vu[1];
            float weightGradXRight = fluid_attributes.dxRights[i].vu[1];
            float weightGradYBottom = fluid_attributes.dyBottoms[i].vu[1];

            float dudx = 0.f;
            float dudy = 0.f;

            dudx += validTopLeftWeight * gridVel_topLeft * (weightGradYTop * fluid_attributes.invSpacing - fluid_attributes.invSpacing);
            dudy += validTopLeftWeight * gridVel_topLeft * (weightGradXLeft * fluid_attributes.invSpacing - fluid_attributes.invSpacing);

            dudx += validTopRightWeight * gridVel_topRight * -(weightGradYTop * fluid_attributes.invSpacing + fluid_attributes.invSpacing);
            dudy += validTopRightWeight * gridVel_topRight * -(weightGradXRight * fluid_attributes.invSpacing);

            dudx += validBottomLeftWeight * gridVel_bottomLeft * -(weightGradYBottom * fluid_attributes.invSpacing);
            dudy += validBottomLeftWeight * gridVel_bottomLeft * -(weightGradXLeft * fluid_attributes.invSpacing + fluid_attributes.invSpacing);

            dudx += validBottomRightWeight * gridVel_bottomRight * (weightGradYBottom * fluid_attributes.invSpacing);
            dudy += validBottomRightWeight * gridVel_bottomRight * (weightGradXRight * fluid_attributes.invSpacing);
        
            fluid_attributes.affineMats[matIdx] = dudx;
            fluid_attributes.affineMats[matIdx + 1] = dudy;
        }
    }

    void transferVGridToParticles(int32_t startIndex, int32_t endIndex) {
        for (int32_t i = startIndex; i < endIndex; ++i) {
            int32_t vi = 2 * i + 1;
            float pvy = fluid_attributes.velocities[vi];
            int matIdx = 4 * i;

            int32_t topLeftCell = topLeftCells[vi];
            int32_t topRightCell = topRightCells[vi];
            int32_t bottomRightCell = bottomRightCells[vi];
            int32_t bottomLeftCell = bottomLeftCells[vi];

            float topLeftWeight = topLeftWeights[vi];
            float topRightWeight = topRightWeights[vi];
            float bottomLeftWeight = bottomLeftWeights[vi];
            float bottomRightWeight = bottomRightWeights[vi];

            float gridVel_topLeft = fluid_attributes.v[topLeftCell];
            float gridVel_topRight = fluid_attributes.v[topRightCell];
            float gridVel_bottomRight = fluid_attributes.v[bottomRightCell];
            float gridVel_bottomLeft = fluid_attributes.v[bottomLeftCell];
           
            const float validTopLeftWeight = fluid_attributes.cellType[topLeftCell] != AIR_CELL || 
                                             fluid_attributes.cellType[topLeftCell - 1] != AIR_CELL;
            const float validTopRightWeight = fluid_attributes.cellType[topRightCell] != AIR_CELL || 
                                              fluid_attributes.cellType[topRightCell - 1] != AIR_CELL;
            const float validBottomLeftWeight = fluid_attributes.cellType[bottomLeftCell] != AIR_CELL || 
                                                fluid_attributes.cellType[bottomLeftCell - 1] != AIR_CELL;
            const float validBottomRightWeight = fluid_attributes.cellType[bottomRightCell] != AIR_CELL || 
                                                 fluid_attributes.cellType[bottomRightCell - 1] != AIR_CELL;

            const float divY = validTopLeftWeight * topLeftWeight + validTopRightWeight * topRightWeight + validBottomLeftWeight * bottomLeftWeight + validBottomRightWeight * bottomRightWeight;

            float picV;
            float corr;
            float flipV;

            if (divY > 0.f) {
                picV = (validTopLeftWeight * topLeftWeight * gridVel_topLeft + 
                        validTopRightWeight * topRightWeight * gridVel_topRight + 
                        validBottomLeftWeight * bottomLeftWeight * gridVel_bottomLeft + 
                        validBottomRightWeight * bottomRightWeight * gridVel_bottomRight
                        ) / divY;

                corr = (validTopLeftWeight * topLeftWeight * (gridVel_topLeft - fluid_attributes.prevV[topLeftCell]) + 
                        validTopRightWeight * topRightWeight * (gridVel_topRight - fluid_attributes.prevV[topRightCell]) + 
                        validBottomLeftWeight * bottomLeftWeight * (gridVel_bottomLeft - fluid_attributes.prevV[bottomLeftCell]) + 
                        validBottomRightWeight * bottomRightWeight * (gridVel_bottomRight - fluid_attributes.prevV[bottomRightCell])
                        ) / divY;

                flipV = pvy + corr;
                fluid_attributes.velocities[vi] = (1.f - fluid_attributes.flipRatio) * picV + fluid_attributes.flipRatio * flipV;
            }

            float weightGradXLeft = fluid_attributes.dxLefts[i].vu[0];
            float weightGradYTop = fluid_attributes.dyTops[i].vu[0];
            float weightGradXRight = fluid_attributes.dxRights[i].vu[0];
            float weightGradYBottom = fluid_attributes.dyBottoms[i].vu[0];

            float dvdx = 0.f;
            float dvdy = 0.f;

            dvdx += validTopLeftWeight * gridVel_topLeft * (weightGradYTop * fluid_attributes.invSpacing - fluid_attributes.invSpacing);
            dvdy += validTopLeftWeight * gridVel_topLeft * (weightGradXLeft * fluid_attributes.invSpacing - fluid_attributes.invSpacing);

            dvdx += validTopRightWeight * gridVel_topRight * -(weightGradYTop * fluid_attributes.invSpacing + fluid_attributes.invSpacing);
            dvdy += validTopRightWeight * gridVel_topRight * -(weightGradXRight * fluid_attributes.invSpacing);

            dvdx += validBottomLeftWeight * gridVel_bottomLeft * -(weightGradYBottom * fluid_attributes.invSpacing);
            dvdy += validBottomLeftWeight * gridVel_bottomLeft * -(weightGradXLeft * fluid_attributes.invSpacing + fluid_attributes.invSpacing);

            dvdx += validBottomRightWeight * gridVel_bottomRight * (weightGradYBottom * fluid_attributes.invSpacing);
            dvdy += validBottomRightWeight * gridVel_bottomRight * (weightGradXRight * fluid_attributes.invSpacing);
        
            fluid_attributes.affineMats[matIdx + 2] = dvdx;
            fluid_attributes.affineMats[matIdx + 3] = dvdy;
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
            int x1 = std::min(x0 + 1, fluid_attributes.numX - 1);

            int y0 = std::max(1, std::min((int)(std::floor((y - halfSpacing) * invSpacing)), fluid_attributes.numY - 2));
            int y1 = std::min(y0 + 1, fluid_attributes.numY - 2);
           
            float dxLeft = ((x - halfSpacing) - x0 * fluid_attributes.cellSpacing) * invSpacing;
            float dyTop = ((y - halfSpacing) - y0 * fluid_attributes.cellSpacing) * invSpacing;
            float dxRight = 1.f - dxLeft;
            float dyBottom = 1.f - dyTop;

            int topLeft = x0 * n + y0;
            int topRight = x1 * n + y0;
            int bottomLeft = x0 * n + y1;
            int bottomRight = x1 * n + y1;

            threadGridBuffers[threadID][topLeft] += dxRight * dyBottom;
            
            threadGridBuffers[threadID][topRight] += dxLeft * dyBottom;
                
            threadGridBuffers[threadID][bottomLeft] += dxRight * dyTop;
                
            threadGridBuffers[threadID][bottomRight] += dxLeft * dyTop;
        }
    }

    void applyLocalDensityUpdates(int start, int end) {
        for (int i = start; i < end; ++i) {
            for (int j = 0; j < fluid_attributes.numThreads; ++j) {
                fluid_attributes.cellDensities[i] += threadGridBuffers[j][i];
            }
        }
    }

    void applyLocalDensityUpdatesMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.gridSize, [this](int start, int end) {
            applyLocalDensityUpdates(start, end);
        });
    }

};