
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

    // Gauss seidel multithreading just isnt viable
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
    }
    
    void drawUVGrids(sf::RenderWindow& window) {
        sf::VertexArray line(sf::Lines, 2);
        int32_t n = numY;
        for (int i = 0; i < 1.f * (WIDTH - 2 * cellSpacing) / cellSpacing; ++i) {
            for (int j = 0; j < 1.f * (HEIGHT - 2 * cellSpacing) / cellSpacing; ++j) {

                // draw u lines (left right)
                float uX = cellSpacing + i * cellSpacing;
                float uY = 1.5 * cellSpacing + j * cellSpacing;
                line[0].position = sf::Vector2f(uX, uY);
                line[0].color  = sf::Color(0, 150, 255);
                line[1].position = sf::Vector2f(uX + u[i * n + j], uY);
                line[1].color = sf::Color(0, 150, 255);
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
    
    void solveIncompressibilityCG() {
        const int32_t n = numY;

        // use doubles for everything when using conjugate gradient; less roundoff error = faster convergence 

        // initial guess at final changes in pressures 
        std::fill(pressure.begin(), pressure.end(), 0.f);

        std::fill(residual.begin(), residual.end(), 0.f);
        std::fill(Ad.begin(), Ad.end(), 0.f);
        std::fill(direction.begin(), direction.end(), 0.f);

        // trying to solve Ap = -b (A * pressure = divergence), where A is the discretized Laplacian, p is the unknown pressure changes in each cell, and b is the initial divergence

        // construct b (set the residual equal to b and the direction equal to residual)
        float scale = 1.f / cellSpacing;
        for (int32_t i = 1; i < numX - 1; ++i) {
            for (int32_t j = 1; j < numY - 1; ++j) {
                const int32_t idx = i * n + j;
                if (cellType[idx] != FLUID_CELL)  {
                    // set non-fluid cells to zero so that their values doesn't mess with the later inner products
                    /*residual[idx] = 0.0f;
                    direction[idx] = 0.0f;
                    Ad[idx] = 0.0f;*/
                    /*continue;
                }

                // residual is -b - Ap, but A starts out as zero because our initial guess at the unknown pressures is just the zero vector. 
                residual[idx] = -scale * (this->u[(i + 1) * n + j] - this->u[idx] + this->v[idx + 1] - this->v[idx]);

                // direction = -b, happens to be equal to residual in this case 
                direction[idx] = residual[idx];
            }
        }

        // iterate
        for (int32_t iter = 0; iter < 1; ++iter) {

            // construct discretized laplacian A_direction (Ad) using the direction vector 

            // we don't actually have to store the gigantic laplacian matrix; think of Ad as the resulting vector of the A * direction matrix multiplication 

            // and direction looks something like this:
                /*
                    [1]
                    [2]
                    [1]
                    [4]
                    [3]
                    [5]
                    [.]
                    [.]
                */

            // The A * d multiplication results in the Ad vector below

            //--------------------------------------------------------------------------------------
            // Extra info:
                // really, A might look something like this (this would be for a 4x4 grid):
                /*  [-4  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0] 
                    [ 1 -4  1  0  0  1  0  0  0  0  0  0  0  0  0  0] 
                    [ 0  1 -4  1  0  0  1  0  0  0  0  0  0  0  0  0] 
                    [ 0  0  1 -4  0  0  0  1  0  0  0  0  0  0  0  0] 
                    [ 1  0  0  0 -4  1  0  0  1  0  0  0  0  0  0  0] 
                    [ 0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0  0] 
                    [ 0  0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0] 
                    [ 0  0  0  1  0  0  1 -4  0  0  0  1  0  0  0  0] 
                    [ 0  0  0  0  1  0  0  0 -4  1  0  0  1  0  0  0] 
                    [ 0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0  0] 
                    [ 0  0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0] 
                    [ 0  0  0  0  0  0  0  1  0  0  1 -4  1  0  0  1] 
                    [ 0  0  0  0  0  0  0  0  1  0  0  1 -4  1  0  0] 
                    [ 0  0  0  0  0  0  0  0  0  1  0  0  1 -4  1  0] 
                    [ 0  0  0  0  0  0  0  0  0  0  1  0  0  1 -4  1] 
                    [ 0  0  0  0  0  0  0  0  0  0  0  1  0  0  1 -4] 

                    The diagonal row of -4s is -1 times the number of non-solid neighbors of cell i (4 on a 2D grid), and every single other cell denoted B_ij is 1 if cell i is a neighbor of cell j, and 0 otherwise

                    The overall structure of this A matrix can only be applied to a 4x4 grid because of the boundary conditions chosen (0 at the sides) 

                        for example, here's A for a 2x2 grid:

                            [-4  1  1  0]
                            [1  -4  0  1]
                            [1   0 -4  1]
                            [0   1  1 -4] 

                        notice how A_4x4[3, 1] != A_16x16[3, 1] because in a 4x4 matrix, those cells are not neighbors, but in a 2x2 matrix, they are

                        The characteristic polynomial of this matrix is (λ+2)(λ+4)^2(λ+6). All values of λ are negative, meaning this matrix is negative definite

                    In this code, a moderate sized grid would be 70x134. The entire laplacian matrix of this grid, A, would require 335.5 MB of memory to store. This is why you shouldn't/can't just solve the matrix using Gaussean elimination/inversion. On top of that, for large, sparse matrices, good iterative methods are faster than direct solvers, not even taking into account cache efficiency. 
                    
                    This matrix A is diagonalizable and only has negative eigenvalues, meaning it's negative-definite. 

                        CG only works for positive-definite matrices, so we just flip the sign of A to make it positive definite, and then swap the signs of A and b, which is why the sign of b is negative (Ap = b ---> -Ap = b ---> Ap = -b)

                    This matrix is also symmetrical, another requirement of CG, meaning that A[i, j] = A[j, i]

                    Why don't we just do p = A^-1 * -b?
                        A is sparse, but its inverse is not. A^-1 would have to be explicitly stored. 
                */
            //--------------------------------------------------------------------------------------

            /*const float invSqrCellSpacing = 1.f / (cellSpacing * cellSpacing);

            for (int32_t i = 1; i < numX - 1; ++i) {
                for (int32_t j = 1; j < numY - 1; ++j) {
                    int32_t idx = i * numY + j;

                    if (cellType[idx] != FLUID_CELL) continue;

                    int32_t leftType = cellType[(i - 1) * n + j] <= AIR_CELL ? 1 : 0;
                    int32_t rightType = cellType[(i + 1) * n + j] <= AIR_CELL ? 1 : 0;
                    int32_t topType = cellType[i * n + j - 1] <= AIR_CELL ? 1 : 0;
                    int32_t bottomType = cellType[i * n + j + 1] <= AIR_CELL ? 1 : 0;

                    // Use Dirichlet boundary conditions for the edges (directly assign them a value of whatever I want; all 0 in this case) 
                    // If a neighbor cell is an air_cell, set its pressure to zero
                    // If a neighbot cell is a solid cell, subtract 1 from the multiplication of the center cell (done automatically with leftType + rightType + topType + bottomType)
                        // also handle solid cells explicitly, which is not what I'm doing 
                    float center = -(leftType + rightType + topType + bottomType) * direction[idx];
                    float left = leftType * direction[(i - 1) * numY + j];
                    float right = rightType * direction[(i + 1) * numY + j];
                    float up = topType * direction[i * numY + j - 1];
                    float down = bottomType * direction[i * numY + j + 1];

                    Ad[idx] = -(center * invSqrCellSpacing + left + right + down + up) * invSqrCellSpacing;
                }
            }

            // α = (residual * residual) / (A_direction * direction)
            // residual * residual
            float alphaNumerator = std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.f);

            // A_direction * direction
            float alphaDenominator = std::inner_product(direction.begin(), direction.end(), Ad.begin(), 0.f);

            float alpha = alphaNumerator / alphaDenominator;

            for (int i = 0; i < numX * numY; ++i) {
                if (cellType[i] != FLUID_CELL) continue;
                pressure[i] += alpha * direction[i];
            }

            for (int i = 0; i < numX * numY; ++i) {
                if (cellType[i] != FLUID_CELL) continue;
                residual[i] -= alpha * Ad[i];
            }

            // β = (residual * residual) / (old_residual * old_residual) 
            float beta = std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.0f) / alphaNumerator;

            for (int i = 0; i < numX * numY; ++i) {
                if (cellType[i] != FLUID_CELL) continue;
                direction[i] = residual[i] + beta * direction[i];
            }
        }

        //std::cout << "\n\n\n";

        // u^(n+1)_(i+1/2, j) = u^(n)_(i+1/2, j) - Δt(p_(i+1, j) - p(i, j))
        // v^(n+1)_(i, j+1/2) = v^(n)_(i, j+1/2) - Δt(p_(i, j+1) - p(i, j))
        scale = dt / cellSpacing;
        for (int i = 1; i < numX - 1; ++i) {
            for (int j = 1; j < numY - 1; ++j) {
                int idx = i * n + j;
                if (cellType[idx] != FLUID_CELL) continue;

                //u[idx] -= scale * (pressure[(i + 1) * n + j] - pressure[idx]);
                //v[idx] -= scale * (pressure[idx + 1] - pressure[idx]);
                u[idx] -= scale * (pressure[idx]);
                v[idx] -= scale * (pressure[idx]);
                u[(i + 1) * n + j] += scale * (pressure[idx]);
                v[idx + 1] += scale * (pressure[idx]);
            }
        }
    }*/
