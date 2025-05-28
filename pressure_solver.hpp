#pragma once
#include "fluid_state.hpp"

struct PressureSolver {
    FluidState &fluid_attributes;
    int n;
    int gridSize;
    int numPressureIters;
    int FLUID_CELL = fluid_attributes.FLUID_CELL;
    int AIR_CELL = fluid_attributes.AIR_CELL;
    int SOLID_CELL = fluid_attributes.SOLID_CELL;

    float k = 8.f;
    float overRelaxation = 1.9f;

    std::vector<double> residual;
    std::vector<double> Adiag;
    std::vector<double> si;
    std::vector<double> li;
    std::vector<double> precon;
    std::vector<double> search;
    std::vector<double> z;
    std::vector<double> dotProducts;
    std::vector<double> pressure;


    PressureSolver(FluidState &fas, int numPressureIters_): fluid_attributes(fas), numPressureIters(numPressureIters_) {
        n = fluid_attributes.numY;

        this->dotProducts.resize(fluid_attributes.numThreads);

        gridSize = fluid_attributes.numX * fluid_attributes.numY;
        this->Adiag.resize(gridSize);
        this->si.resize(gridSize);
        this->li.resize(gridSize);
        this->precon.resize(gridSize);
        this->z.resize(gridSize);
        this->search.resize(gridSize);
        this->residual.resize(gridSize);
        this->pressure.resize(gridSize);
    }

    void setUpResidual() {
        std::fill(begin(residual), end(residual), 0.0);

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

    void ScaledAdd(std::vector<double>& a, std::vector<double>& b, double c) {
        for (int i = 0; i < gridSize; ++i) {
            if (fluid_attributes.cellType[i] == FLUID_CELL) {
                a[i] += b[i] * c;
            }
        }
    }

    double DotMulti(std::vector<double> *a, std::vector<double> *b) {
        const int32_t numThreads_ = 1;
        const int32_t numColumnsPerThread = (gridSize) / numThreads_;
        const int32_t numMissedColumns = gridSize - numColumnsPerThread * numThreads_;

        std::fill(begin(dotProducts), end(dotProducts), 0.0);

        for (int i = 0; i < numThreads_; ++i) {
            int start = i * numColumnsPerThread;
            int end = (i == numThreads_ - 1) ? (gridSize) : (start + numColumnsPerThread);
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
        const int32_t numColumnsPerThread = gridSize / numThreads_;
        const int32_t numMissedColumns = gridSize - numColumnsPerThread * numThreads_;

        for (int i = 0; i < numThreads_; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->EqualsPlusTimes(a, b, c, i * numColumnsPerThread, i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->EqualsPlusTimes(a, b, c, gridSize - numMissedColumns, gridSize);

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

    void setUpA() {
        float scale = 1;//dt / (cellSpacing * cellSpacing); // also see what happens when you divide by 1000 (density of water) as well

        std::fill(Adiag.begin(), Adiag.end(), 0.0);

        // Ax[i] is the cell to the right, Ay[i] is the cell below
        std::fill(si.begin(), si.end(), 0.0);
        std::fill(li.begin(), li.end(), 0.0);

        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int idx = i * n + j;

                if (fluid_attributes.cellType[idx] != fluid_attributes.FLUID_CELL) continue;

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

        for (int i = 1; i < fluid_attributes.numX - 1; ++i) { // I call this the level iteration. Every time j increases, we are at the next level iteration
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) { // and this is the swipe iteration. Every time i increases, we are at the next swipe iteration
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
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
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

        for (int i = fluid_attributes.numX - 2; i > 0; --i) {
            for (int j = fluid_attributes.numY - 2; j > 0; --j) {
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
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
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

    void projectPCG(int numIters) {
        std::fill(begin(pressure), end(pressure), 0.f);

        setUpResidual();
        setUpA();
        buildPreconditioner();
        applyPreconditioner(&z, &residual);
        std::copy(begin(z), end(z), begin(search));

        // DotMulti needs to be debugged, EqualsPlusTimesMulti good

        /*double sigma = 0.0;
        Dot(&z, &residual, 0, numX * numY, sigma);*/

        double sigma = DotMulti(&z, &residual);

        for (int iter = 0; iter < numIters && sigma > 0; ++iter) {
            matVec(&z, &search);

            //double denom = 0.0;
            //Dot(&z, &search, 0, numX * numY, denom);
            //double alpha = sigma / denom;

            double alpha = sigma / DotMulti(&z, &search);

            ScaledAdd(pressure, search, alpha);
            ScaledAdd(residual, z, -alpha);

            applyPreconditioner(&z, &residual);  // applying the preconditioner is the only thing making MICCG(0) so slow

            //double sigmaNew = 0.0;
            //Dot(&z, &residual, 0, numX * numY, sigmaNew);

            double sigmaNew = DotMulti(&z, &residual);
            
            EqualsPlusTimesMulti(&search, &z, sigmaNew / sigma);

            sigma = sigmaNew;
        }

        applyPressureMulti();
    }

    void applyPressureMulti() {
        const int32_t numThreads_ = 1;
        const int32_t numColumnsPerThread = (fluid_attributes.numX - 2) / numThreads_;
        const int32_t numMissedColumns = (fluid_attributes.numX - 2) - numColumnsPerThread * numThreads_;

        for (int i = 0; i < numThreads_; ++i) {
            fluid_attributes.thread_pool.addTask([&, i]() {
                this->applyPressure(1 + i * numColumnsPerThread, 1 + i * numColumnsPerThread + numColumnsPerThread);
            });
        }

        this->applyPressure(fluid_attributes.numX - 1 - numMissedColumns, fluid_attributes.numX - 1);

        fluid_attributes.thread_pool.waitForCompletion();

        /*const int32_t numColumns = (fluid_attributes.numX - 2);
        fluid_attributes.thread_pool.dispatch(numColumns, [this](int start, int end) {
            this->applyPressure(start, end);
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

    void projectSOR(int numIters) {
        for (int iter = 0; iter < numIters; ++iter) {
            for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
                for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
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
    }

    void passRedBlackGS(int start, int stop, bool red) {
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

    void projectRedBlackGS(int numIters) {
        for (int i = 0; i < numIters; ++i) {
            passRedBlackGS(0, fluid_attributes.numX - 1, 0);
            passRedBlackGS(0, fluid_attributes.numX - 1, 1);
        }
    }

    void projectRedBlackGSMulti(int numIters, int numThreads) {
        int columnsPerThread = (fluid_attributes.numX - 1) / numThreads;
        int numMissedColumns = fluid_attributes.numX - 1 - numThreads * columnsPerThread;

        // not 100% thread safe but who cares there's no artifacts and it runs much faster than thread safe version
        for (int iter = 0; iter < numIters; ++iter) {
            for (int i = 0; i < numThreads; ++i) {
                fluid_attributes.thread_pool.addTask([this, columnsPerThread, i] {
                    int start = i * columnsPerThread + 1;
                    int end = i * columnsPerThread + columnsPerThread + 1;
                
                    passRedBlackGS(start, end, 0);
                    passRedBlackGS(start, end, 1);
                });
            }
        
            passRedBlackGS(fluid_attributes.numX - 1 - numMissedColumns, fluid_attributes.numX - 1, 0);
            passRedBlackGS(fluid_attributes.numX - 1 - numMissedColumns, fluid_attributes.numX - 1, 1);
        
            fluid_attributes.thread_pool.waitForCompletion();
        }
    }

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
