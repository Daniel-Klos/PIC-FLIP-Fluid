    // MICCG(0) code, as described in Bridson's book
    // --------------------------------------------------------------------------------------------------------------------------
    void setUpA() {
        float scale = 1;//dt / (density * cellSpacing * cellSpacing);

        std::fill(Adiag.begin(), Adiag.end(), 0.0);

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

        for (int i = 1; i < fluid_attributes.numX - 1; ++i) { // I call this the level iteration. Every time i increases, we are at the next level iteration
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) { // and this is the swipe iteration. Every time j increases, we are at the next swipe iteration
                int idx = i * n + j;

                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;

                double e = Adiag[idx];

                // si is short for swipe iteration, just some slang for you new gens
                    // si[index] = 0 if the cell in the next swipe iteration touching the current is not fluid, else si[index] = -1

                // li is short for level iteration
                    // li[index] = 0 if the cell in the next level iteration touching the current cell is not fluid, else li[index] = -1

                // if I add/subtract something by S in a comment, that means I'm taking the value of it at the next/previous swipe iteration (exactly one row below/above, one column left/right)

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

                // don't encase a small amount of fluid cells within solids, it messes with the preconditioner. But modifying it only if e is big enough patches up the problem
                if (e * e > 1e-18) {
                    precon[idx] = 1.0 / std::sqrt(e);
                }
            }
        }
    }

    void applyPreconditioner(std::vector<double> &dst, std::vector<double> &a) {
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int32_t idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;

                double t = a[idx];

                // if cell - S is fluid, then do: (si - S) * (precon - S) * (dst - S)
                if (fluid_attributes.cellType[idx - 1] == FLUID_CELL) {
                    t -= si[idx - 1] * precon[idx - 1] * dst[idx - 1];
                }

                // if cell - L is fluid, then do: (si - L) * (precon - L) * (dst - L)
                if (fluid_attributes.cellType[idx - n] == FLUID_CELL) {
                    t -= li[idx - n] * precon[idx - n] * dst[idx - n];
                }

                dst[idx] = t * precon[idx];
            }
        }

        for (int i = fluid_attributes.numX - 2; i > 0; --i) {
            for (int j = fluid_attributes.numY - 2; j > 0; --j) {
                int32_t idx = i * n + j;
                if (fluid_attributes.cellType[idx] != FLUID_CELL) continue;

                double t = dst[idx];

                // if cell + S is fluid, then do: (si) * (precon) * (dst + S)
                if (fluid_attributes.cellType[idx + 1] == FLUID_CELL) {
                    t -= si[idx] * precon[idx] * dst[idx + 1];
                }

                // if cell + L is fluid, then do: (li) * (precon) * (dst + L)
                if (fluid_attributes.cellType[idx + n] == FLUID_CELL) {
                    t -= li[idx] * precon[idx] * dst[idx + n];
                }

                dst[idx] = t * precon[idx];
            }
        }
    }

    void matVec(std::vector<double> &dst, std::vector<double> &b) {
        for (int i = 1; i < fluid_attributes.numX - 1; ++i) {
            for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
                int32_t idx = i * n + j;
                
                double t = Adiag[idx] * b[idx];

                // (si - S) * (b - S)
                t += si[idx - 1] * b[idx - 1];

                // (li - L) * (b - L)
                t += li[idx - n] * b[idx - n];

                // (si) * (b + S)
                t += si[idx] * b[idx + 1];

                // (li) * (b + L)
                t += li[idx] * b[idx + n];

                dst[idx] = t;
            }
        }
    }

    void projectMICCG(int numIters) {
        std::fill(begin(pressure), end(pressure), 0.f);

        setUpResidual();
        setUpA();
        buildPreconditioner();
        applyPreconditioner(z, residual);
        std::copy(begin(z), end(z), begin(direction));

        double sigma = DotMulti(z, residual);

        for (int iter = 0; iter < numIters && sigma > 0; ++iter) {
            matVec(z, direction);

            double alpha = sigma / DotMulti(z, direction);

            ScaledAdd(pressure, direction, alpha);
            ScaledAdd(residual, z, -alpha);

            applyPreconditioner(z, residual);

            double sigmaNew = DotMulti(z, residual);
            
            EqualsPlusTimesMulti(direction, z, sigmaNew / sigma);

            sigma = sigmaNew;
        }

        applyPressureMulti();
    }
    // --------------------------------------------------------------------------------------------------------------------------

    // constant memory collision detection/resolution code
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
    
    void drawHashNums(sf::RenderWindow& window, sf::Text& text) {
        for (int i = 0; i < WIDTH; i += this->spacing) {
            for (int j = 0; j < HEIGHT; j += this->spacing) {
                text.setPosition(i, j);
                text.setString(std::to_string(this->hashCoords(i, j)));
                window.draw(text);
            }
        }
    }
    
    
    int hashCoords(int xi, int yi) {
        int hash = (intCoord(xi) * 92837111) ^ (intCoord(yi) * 689287499);
        return std::abs(hash) % this->tableSize;
    }

    int intCoord(int coord) {
        return std::floor(coord / this->spacing) + 1;
    }


    void initializePhiGrid() {
        const float positive = 0;
        for (int i = 0; i < numX; ++i) {
            for (int j = 0; j < numY; ++j) {
                if (cellType[i * numY + j] == SOLID_CELL) {
                    phiGrid[i * numY + j] = -1.f;
                }
            }
        }
    }
    
    void fastSweepLoop() {
        for (int i = 0; i < numX; ++i) {
            for (int j = 0; j < numY; ++j) {
                if (cellType[i * n + j] != SOLID_CELL) continue; 
    
                float phiX2 = std::max(phiGrid[std::max(i - 1, 0) * numY + j], phiGrid[std::min(i + 1, numX - 1) * numY + j]);
                float phiY2 = std::max(phiGrid[i * numY + std::max(j - 1, 0)], phiGrid[i * numY + std::min(j + 1, numY - 1)]);
    
                float newPhi = std::min(phiGrid[i * numY + j], -1.0f + std::max(phiX2, phiY2));
                phiGrid[i * numY + j] = newPhi;
            }
        }
    }
    
    void fastSweep() {
        fastSweepLoop();
        fastSweepLoop();
        fastSweepLoop();
        fastSweepLoop();

        /*fastSweepLoop(0, numX, 0, numY, 1, 1);
        fastSweepLoop(numX - 1, -1, numY - 1, -1, -1, -1);
        fastSweepLoop(0, numX, numY - 1, -1, 1, -1);
        fastSweepLoop(numX - 1, -1, 0, numY, -1, 1);*/
    }

    void computeSDF() {
        initializePhiGrid();
        fastSweep();
    }

    float interpolatePhi(float px, float py) {
        float i = std::floor((px - cellSpacing / 2) * invSpacing);
        float j = std::floor((py - cellSpacing / 2) * invSpacing);

        int i0 = std::clamp(static_cast<int>(i), 0, numX - 1);
        int j0 = std::clamp(static_cast<int>(j), 0, numY - 1);
        int i1 = std::clamp(i0 + 1, 0, numX - 1);
        int j1 = std::clamp(j0 + 1, 0, numY - 1);

        float fx = (px - cellSpacing / 2.f - cellSpacing * i) * invSpacing;
        float fy = (py - cellSpacing / 2.f - cellSpacing * j) * invSpacing;
        float sx = 1.f - fx;
        float sy = 1.f - fy;

        i = static_cast<int>(i);
        j = static_cast<int>(j);

        float phi00 = phiGrid[i0 * n + j0];  // top left
        float phi10 = phiGrid[i1 * n + j0];  // top right
        float phi01 = phiGrid[i0 * n + j1];  // bottom left
        float phi11 = phiGrid[i1 * n + j1];  // bottom right

        return sx * sy * phi00 + 
               fx * sy * phi10 + 
               fy * sx * phi01 + 
               fx * fy * phi11;
    }

    void computeNormal(float &nx, float& ny, float px, float py, float phi) {
        // make it so that it only computes diagonal normals when interpolatePhi() > 1 or something
        int i = static_cast<int>(px / cellSpacing);
        int j = static_cast<int>(py / cellSpacing);

        int left = std::max(i - 1, 0);
        int right = std::min(i + 1, numX - 1);
        int top = std::max(j - 1, 0);
        int bottom = std::min(j + 1, numY - 1);

        float dPhidx = (phiGrid[right * n + j] - phiGrid[left * n + j]) / (2.f * cellSpacing);
        float dPhidy = (phiGrid[i * n + bottom] - phiGrid[i * n + top]) / (2.f * cellSpacing);

        
        float dx = px - i * cellSpacing;
        float dy = py - j * cellSpacing;

        float length = std::sqrt(dPhidx * dPhidx + dPhidy * dPhidy);
        if (length > 0 && length < WIDTH) {
            nx = dPhidx / length;
            ny = dPhidy / length;
        }
        else {
            nx = ny = 0;
        }
    }

    void separateParticle(float& px, float& py, float& vx, float& vy) {
        float phi = interpolatePhi(px, py);
        if (phi < 0) {
            float nx, ny;
            computeNormal(nx, ny, px, py, phi);
            px += -phi * nx * 10;
            py += -phi * ny * 10;

            float velocityNormal = vx * nx + vy * ny;
            if (velocityNormal < 0) {
                vx -= velocityNormal * nx;
                vy -= velocityNormal * ny;
            }
        }
    }

    void collideSurfaces(const uint32_t start, const uint32_t end) {
        for (int i = start; i < end; ++i) {
            separateParticle(positions[2 * i], positions[2 * i + 1], velocities[2 * i], velocities[2 * i + 1]);
        }
    }

    void showSeparationMouse(sf::RenderWindow& window) {
        float phi = interpolatePhi(mouseX, mouseY);
        if (phi < 0) {
            float nx, ny;
            computeNormal(nx, ny, mouseX, mouseY, phi);

            float pointX = mouseX - phi * nx;
            float pointY = mouseY - phi * ny;

            sf::VertexArray line(sf::Lines, 2);
            line[0].position = sf::Vector2f(mouseX, mouseY);
            line[0].color  = sf::Color(255, 0, 0);
            line[1].position = sf::Vector2f(mouseX - phi * nx * 100, mouseY - phi * ny * 100);
            line[1].color = sf::Color(255, 0, 0);
            window.draw(line);
        }
    }

    void drawPhiValues(sf::RenderWindow& window) {
        float px = cellSpacing / 3;
        float py = cellSpacing / 3;
        for (int i = 0; i < numX; ++i) {
            for (int j = 0; j < numY; ++j) {
                text.setString(std::to_string(static_cast<int>(phiGrid[i * n + j])));
                text.setPosition(px, py);
                window.draw(text);
                py += cellSpacing;
            }
            py = cellSpacing / 3;
            px += cellSpacing;
        }
    }
