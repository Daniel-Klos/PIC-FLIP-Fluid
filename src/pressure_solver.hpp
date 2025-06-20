#pragma once
#include <functional>
#include "simulation_state.hpp"

struct PressureSolver {
    FluidState &fluid_attributes;
    int n;
    int gridSize;
    int numPressureIters = 30;
    int FLUID_CELL = fluid_attributes.FLUID_CELL;
    int AIR_CELL = fluid_attributes.AIR_CELL;
    int SOLID_CELL = fluid_attributes.SOLID_CELL;

    float k = 10.f;
    float overRelaxation = 1.9f;

    std::vector<double> residual;
    std::vector<double> Adiag;
    std::vector<double> precon;
    std::vector<double> direction;
    std::vector<double> dotProducts;
    std::vector<double> pressure;

    std::vector<uint8_t> neighbors;
    static constexpr uint8_t LEFT = 64;
    static constexpr uint8_t RIGHT = 32;
    static constexpr uint8_t BOTTOM = 16;
    static constexpr uint8_t TOP = 8;
    static constexpr uint8_t CENTER = 7;
    static constexpr std::array<uint8_t, 4> directions{LEFT, RIGHT, BOTTOM, TOP};

    std::vector<float> coarseU;
    int max_multigrid_coarseness = 4;

    PressureSolver(FluidState &fas): fluid_attributes(fas) {
        n = fluid_attributes.numY;

        this->dotProducts.resize(fluid_attributes.numThreads);

        gridSize = fluid_attributes.numX * fluid_attributes.numY;
        this->Adiag.resize(gridSize);
        this->neighbors.resize(gridSize);
        this->precon.resize(gridSize);
        this->direction.resize(gridSize);
        this->residual.resize(gridSize);
        this->pressure.resize(gridSize);
    }


    // What I'm using currently
    // --------------------------------------------------------------------------------------------------------------------------------------------

    void passRedBlackSOR(int start, int stop, bool red) {
        for (int i = start; i < stop; ++i) {
            for (int j = (i + red) % 2; j < fluid_attributes.numY - 1; j += 2) {
                int idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;

                float leftType = fluid_attributes.cellType[idx - n] <= AIR_CELL ? 1 : 0;
                float rightType = fluid_attributes.cellType[idx + n] <= AIR_CELL ? 1 : 0;
                float topType = fluid_attributes.cellType[idx - 1] <= AIR_CELL ? 1 : 0;
                float bottomType = fluid_attributes.cellType[idx + 1] <= AIR_CELL ? 1 : 0;

                float divideBy = leftType + rightType + topType + bottomType;
            
                if (divideBy == 0.f) continue;

                float divergence = fluid_attributes.u[idx + n] - fluid_attributes.u[idx] + fluid_attributes.v[idx + 1] - fluid_attributes.v[idx];
                if (fluid_attributes.particleRestDensity > 0.f) {
                    float compression = fluid_attributes.cellDensities[idx] - fluid_attributes.particleRestDensity;
                    if (compression > 0.f) {
                        divergence -= k * compression;
                    }
                }

                float p = divergence / divideBy;
                p *= overRelaxation;

                fluid_attributes.u[idx] += leftType * p;
                fluid_attributes.u[idx + n] -= rightType * p;
                fluid_attributes.v[idx] += topType * p;
                fluid_attributes.v[idx + 1] -= bottomType * p;
            }
        }
    }

    void projectRedBlackSORMulti(int numIters) {
        int columnsPerThread = (fluid_attributes.numX - 1) / fluid_attributes.numThreads;
        int numMissedColumns = fluid_attributes.numX - 1 - fluid_attributes.numThreads * columnsPerThread;

        for (int iter = 0; iter < numIters; ++iter) {
            for (int i = 0; i < fluid_attributes.numThreads; ++i) {
                fluid_attributes.thread_pool.addTask([this, columnsPerThread, i] {
                    int start = i * columnsPerThread + 1;
                    int end = i * columnsPerThread + columnsPerThread + 1;
                
                    passRedBlackSOR(start, end, 0);
                });
            }
        
            passRedBlackSOR(fluid_attributes.numX - numMissedColumns, fluid_attributes.numX - 1, 0);
        
            fluid_attributes.thread_pool.waitForCompletion();

            for (int i = 0; i < fluid_attributes.numThreads; ++i) {
                fluid_attributes.thread_pool.addTask([this, columnsPerThread, i] {
                    int start = i * columnsPerThread + 1;
                    int end = i * columnsPerThread + columnsPerThread + 1;
                
                    passRedBlackSOR(start, end, 1);
                });
            }
        
            passRedBlackSOR(fluid_attributes.numX - numMissedColumns, fluid_attributes.numX - 1, 1);
        
            fluid_attributes.thread_pool.waitForCompletion();
        }
    }
    // --------------------------------------------------------------------------------------------------------------------------------------------


    // Soon to be Multigrid code
    // --------------------------------------------------------------------------------------------------------------------------------------------
    void restrict() {
        
    }

    void prolongate() {

    }

    void projectVCycle() {

    }
    // --------------------------------------------------------------------------------------------------------------------------------------------



    // conjugate gradient code
    // --------------------------------------------------------------------------------------------------------------------------------------------

    /*void setUpResidualMulti() {
        std::fill(begin(residual), end(residual), 0.0);

        fluid_attributes.thread_pool.dispatch(fluid_attributes.numColumns, [this])
    }*/

    void setUpResidual() {
        for (int32_t i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int32_t j = 1; j < fluid_attributes.numY - 1; ++j) {
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

    void ScaledAddMulti(std::vector<double> &a, std::vector<double> &b, double c) {
        const int32_t numCellsPerThread = fluid_attributes.gridSize / fluid_attributes.numThreads;
        const int32_t numMissedColumns = fluid_attributes.gridSize - numCellsPerThread * fluid_attributes.numThreads;

        fluid_attributes.thread_pool.dispatch(fluid_attributes.gridSize, [this, &a, &b, c](int start, int end) {
            ScaledAdd(a, b, c, start, end);
        });
    }

    void ScaledAdd(std::vector<double> &a, std::vector<double> &b, double c, int start, int end) {
        for (int i = start; i < end; ++i) {
            if (fluid_attributes.cellType[i] == FLUID_CELL) {
                a[i] += b[i] * c;
            }
        }
    }

    double DotMulti(std::vector<double> &a, std::vector<double> &b) {
        const int32_t numCellsPerThread = fluid_attributes.gridSize / fluid_attributes.numThreads;
        const int32_t numMissedColumns = fluid_attributes.gridSize - numCellsPerThread * fluid_attributes.numThreads;

        std::fill(begin(dotProducts), end(dotProducts), 0.0);

        for (int i = 0; i < fluid_attributes.numThreads; ++i) {
            int start = i * numCellsPerThread;
            int end = (i == fluid_attributes.numThreads - 1) ? (gridSize) : (start + numCellsPerThread);
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

    void Dot(std::vector<double> &a, std::vector<double> &b, int start, int end, double &res) {
        for (int i = start; i < end; ++i) {
            if (fluid_attributes.cellType[i] == FLUID_CELL) {
                res += a[i] * b[i];
            }
        }
    }

    void EqualsPlusTimesMulti(std::vector<double> &a, std::vector<double> &b, double c) {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.gridSize, [this, &a, &b, c](int start, int end) {
            EqualsPlusTimes(a, b, c, start, end);
        });
    }

    void EqualsPlusTimes(std::vector<double> &a, std::vector<double> &b, double c, int start, int end) {
        for (int i = start; i < end; ++i) {
            if (fluid_attributes.cellType[i] == FLUID_CELL) {
                a[i] = b[i] + a[i] * c;
            }
        }
    }

    void applyPressureMulti() {
        const int32_t numColumnsPerThread = (fluid_attributes.numX - 2) / fluid_attributes.numThreads;
        const int32_t numMissedColumns = (fluid_attributes.numX - 2) - numColumnsPerThread * fluid_attributes.numThreads;

        for (int i = 0; i < fluid_attributes.numThreads; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->applyPressure(1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->applyPressure(fluid_attributes.numX - 1 - numMissedColumns, fluid_attributes.numX - 1);

        fluid_attributes.thread_pool.waitForCompletion();

        /*const int32_t numColumns = (fluid_attributes.numX - 2);
        fluid_attributes.thread_pool.dispatch(numColumns, [this](int start, int end) {
            fluid_attributes.applyPressure(start, end);
        });*/
    }

    void applyPressure(int start, int end) {
        //const float density = 1000.f;
        //const float scale = dt / (density * cellSpacing);
    
        for (int i = start; i < end; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
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

    uint8_t getMaterialFromDirection(uint8_t direction, int idx) {
        switch (direction) {
            case LEFT:
                return fluid_attributes.cellType[idx - n];
            case RIGHT:
                return fluid_attributes.cellType[idx + n];
            case BOTTOM:
                return fluid_attributes.cellType[idx + 1];
            case TOP:
                return fluid_attributes.cellType[idx - 1];
        }

        return SOLID_CELL; // should never get executed
    }

    uint8_t updateNbrFromNeighbor(uint8_t material, uint8_t nbr_info, uint8_t dir) {

        // if neighbor cell is fluid or air, add to the first 3 bits
        if (material != SOLID_CELL) {
            nbr_info++;
        }

        // if neighbor cell is not fluid, end here
        if (material != FLUID_CELL) {
            return nbr_info;
        }

        // else, add the direction onto the last 4 bits of nbr_info. Using OR here because we don't wanna touch the first 3 bits
        return nbr_info | dir;
    }

    void setUpNeighbors() {
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;

                uint8_t nbr_info = 0u;
                for (uint8_t dir : directions) {
                    uint8_t material = getMaterialFromDirection(dir, idx);
                    nbr_info = updateNbrFromNeighbor(material, nbr_info, dir);
                }

                neighbors[idx] = nbr_info;
            }
        }
    }

    void MatVecMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.numX - 2, [this](int start, int end) {
            MatVec(start, end);
        });
    }

    void MatVec(int start, int end) {
        for (int i = start + 1; i < end + 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int idx = i * n + j;

                uint8_t nbrs = neighbors[idx];
                if (!nbrs) { // if a cell has no neighbors continue
                    Adiag[idx] = 0.0;
                    continue;
                }

                // all of these & statements are just looking at different bits in nbrs
                // nbrs & CENTER just extracts the first 3 bits of nbrs; if nbrs is 1001 010, then nbrs & center = 1001 010 & 0000 111 = 0000 010
                Adiag[idx] = 
                    ((nbrs & CENTER) * direction[idx]) -
                    ((nbrs & LEFT) ? direction[idx - n] : 0) -
                    ((nbrs & RIGHT) ? direction[idx + n] : 0) -
                    ((nbrs & BOTTOM) ? direction[idx + 1] : 0) -
                    ((nbrs & TOP) ? direction[idx - 1] : 0);
            }
        }
    }

    // conjugate gradient
    void projectCG() {

        setUpNeighbors();

        std::fill(begin(pressure), end(pressure), 0.0);
        std::fill(begin(Adiag), end(Adiag), 0.0);

        setUpResidual();

        std::copy(begin(residual), end(residual), begin(direction));

        double sigma = DotMulti(residual, residual);
        
        for (int iter = 0; iter < numPressureIters && sigma > 0; ++iter) {
            MatVecMulti();
            
            double alpha = sigma / DotMulti(direction, Adiag);
            ScaledAddMulti(pressure, direction, alpha);
            ScaledAddMulti(residual, Adiag, -alpha);
            
            double sigmaOld = sigma;
            sigma = DotMulti(residual, residual);
            EqualsPlusTimesMulti(direction, residual, sigma / sigmaOld);
        }

        applyPressureMulti();
    }

    // general preconditioned conjugate gradient (using matrixless GS/Jacobi methods, so no z needed)
    void projectPCG(std::function<void()> precondition) {

        setUpNeighbors();

        std::fill(begin(pressure), end(pressure), 0.0);
        std::fill(begin(Adiag), end(Adiag), 0.0);

        setUpResidual();

        precondition();

        std::copy(begin(residual), end(residual), begin(direction));

        double sigma = DotMulti(residual, residual);
        
        for (int iter = 0; iter < numPressureIters && sigma > 0; ++iter) {
            MatVecMulti();
            
            double alpha = sigma / DotMulti(direction, Adiag);
            ScaledAddMulti(pressure, direction, alpha);
            ScaledAddMulti(residual, Adiag, -alpha);

            precondition();
            
            double sigmaOld = sigma;
            sigma = DotMulti(residual, residual);
            EqualsPlusTimesMulti(direction, residual, sigma / sigmaOld);
        }

        applyPressureMulti();
    }
    // --------------------------------------------------------------------------------------------------------------------------------------------


    void addToNumPressureIters(int32_t add) {
        numPressureIters += add;
    }

    int32_t getNumPressureIters() {
        return numPressureIters;
    }

    void addToDivergenceModifier(float add) {
        this->k += add;
    }

    float getDivergenceModifier() {
        return this->k;
    }
};