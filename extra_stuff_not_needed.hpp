
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

    /*
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

        /*
        float dx = px - i * cellSpacing;
        float dy = py - j * cellSpacing;

        */
/*
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
        //std::cout << phi << "\n";
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
    }*/
