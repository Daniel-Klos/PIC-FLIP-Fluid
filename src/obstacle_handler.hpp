#pragma once
#include "simulation_state.hpp"
#include "obstacle_rendering.hpp"
#include "utils.hpp"

struct ObstacleHandler {
    FluidState &fluid_attributes;
    ObstacleRenderer &obstacle_renderer;

    int pencilRadius = 1;

    bool solidDrawing = false;

    ObstacleHandler(FluidState &fas, ObstacleRenderer &ore): fluid_attributes(fas), obstacle_renderer(ore) {}

    void constrainWalls(const uint32_t startIndex, const uint32_t endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            int xi = 2 * i;
            int yi = xi + 1;
            if (fluid_attributes.positions[xi] - fluid_attributes.radius < fluid_attributes.cellSpacing) {
                fluid_attributes.positions[xi] = fluid_attributes.radius + fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[xi] < 0) {
                    fluid_attributes.velocities[xi] = 0.f;
                }
            }
            else if (fluid_attributes.positions[xi] + fluid_attributes.radius > fluid_attributes.frame_context.WIDTH - fluid_attributes.cellSpacing) {
                fluid_attributes.positions[xi] = fluid_attributes.frame_context.WIDTH - fluid_attributes.radius - fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[xi] > 0) {
                    fluid_attributes.velocities[xi] = 0.f;
                }
            }
            if (fluid_attributes.positions[yi] - fluid_attributes.radius < fluid_attributes.cellSpacing) {
                fluid_attributes.positions[yi] = fluid_attributes.radius + fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[yi] < 0) {
                    fluid_attributes.velocities[yi] = 0.f;
                }
            }
            else if (fluid_attributes.positions[yi] + fluid_attributes.radius > fluid_attributes.frame_context.HEIGHT - fluid_attributes.cellSpacing) {
                fluid_attributes.positions[yi] = fluid_attributes.frame_context.HEIGHT - fluid_attributes.radius - fluid_attributes.cellSpacing;
                if (fluid_attributes.velocities[yi] > 0) {
                    fluid_attributes.velocities[yi] = 0.f;
                }
            }
        }
    }

    void constrainWallsMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->constrainWalls(start, end);
        });
    }

    float calculateBoxNormals(float bx, float by, float px, float py, float& nx, float& ny) {
        float localpx = px - bx;
        float localpy = py - by;
        float dx = std::abs(localpx) - fluid_attributes.halfSpacing;
        float dy = std::abs(localpy) - fluid_attributes.halfSpacing;
        float pdx = std::max(dx, 0.f);
        float pdy = std::max(dy, 0.f);
        float insideDist = std::min(std::max(dx, dy), 0.f);
        float dist = sqrt(pdx * pdx + pdy * pdy) + insideDist - fluid_attributes.radius;
    
        if (dist >= 0.f) {
            return dist;
        }
    
        float dirX = sign(localpx);
        float dirY = sign(localpy);

        int idx = fluid_attributes.n * static_cast<int>(px * fluid_attributes.invSpacing) + static_cast<int>(py * fluid_attributes.invSpacing);
    
        if (pdx > 0 && pdy > 0 && fluid_attributes.cellType[idx + dirX] != fluid_attributes.SOLID_CELL && fluid_attributes.cellType[idx + dirY] != fluid_attributes.SOLID_CELL) {
            float len = dist - (insideDist - fluid_attributes.radius);
            nx = -dirX * dx / len;
            ny = -dirY * dy / len;
        }
        else if (abs(localpx) > abs(localpy)) {
            nx = dirX * sign(dist);
            ny = 0;
        }
        else {
            nx = 0;
            ny = dirY * sign(dist);
        }
    
        return dist;
    }

    void separateParticle(float& px, float& py, float& vx, float& vy) {
        int localX = px / fluid_attributes.cellSpacing;
        int localY = py / fluid_attributes.cellSpacing;

        int x0 = std::max(0, localX - 1);
        int x1 = std::min(fluid_attributes.numX - 1, localX + 1);
        int y0 = std::max(0, localY - 1);
        int y1 = std::min(fluid_attributes.numY - 1, localY + 1);

        float prevX = px;
        float prevY = py;

        for (int i = x0; i <= x1; ++i) {
            for (int j = y0; j <= y1; ++j) {
                if (fluid_attributes.cellType[i * fluid_attributes.n + j] == fluid_attributes.SOLID_CELL) {
                    float nx = 0;
                    float ny = 0;
                    float dist = calculateBoxNormals(i * fluid_attributes.cellSpacing + fluid_attributes.halfSpacing, j * fluid_attributes.cellSpacing + fluid_attributes.halfSpacing, px, py, nx, ny);
                    if (dist < 0.f) {
                        float velocityNormal = vx * -nx + vy * -ny;
                        if (velocityNormal < 0) {
                            vx -= velocityNormal * -nx;
                            vy -= velocityNormal * -ny;
                        }
                        px += nx * dist;
                        py += ny * dist;
                    }
                }
            }
        }
    }

    void collideSurfaces(const uint32_t start, const uint32_t end) {
        for (int i = start; i < end; ++i) {
            separateParticle(fluid_attributes.positions[2 * i], fluid_attributes.positions[2 * i + 1], fluid_attributes.velocities[2 * i], fluid_attributes.velocities[2 * i + 1]);
        }
    }

    void collideSurfacesMulti() {
        fluid_attributes.thread_pool.dispatch(fluid_attributes.num_particles, [this](int start, int end) {
            this->collideSurfaces(start, end);
        });
    }

    void drawSolids() {
        int localX = static_cast<int>(fluid_attributes.frame_context.world_mouse_pos.x / fluid_attributes.cellSpacing);
        int localY = static_cast<int>(fluid_attributes.frame_context.world_mouse_pos.y / fluid_attributes.cellSpacing);

        int numObstacles = fluid_attributes.obstaclePositions.size();

        int numPotentialAddedObstacles = ((2 * pencilRadius + 1) * (2 * pencilRadius + 1));
        fluid_attributes.obstaclePositions.resize(numObstacles + numPotentialAddedObstacles);

        int numAddedObstacles = 0;
        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                int x = localX + i;
                int y = localY + j;

                if (x > 0 && y > 0 && x < fluid_attributes.numX - 1 && y < fluid_attributes.numY - 1) {

                    int idx = x * fluid_attributes.numY + y;

                    if (fluid_attributes.cellType[idx] != fluid_attributes.SOLID_CELL) {
                        fluid_attributes.cellType[idx] = fluid_attributes.SOLID_CELL;
                        fluid_attributes.obstaclePositions[numObstacles + numAddedObstacles] = sf::Vector2i{x, y};
                        numAddedObstacles++;
                    }
                }
            }
        }

        fluid_attributes.obstaclePositions.resize(numObstacles + numAddedObstacles);

        if (numAddedObstacles > 0) {
            obstacle_renderer.obstacleVa.resize(4 * (numObstacles + numAddedObstacles));

            sf::Color gray = sf::Color(150, 150, 150);
            for (int i = 0; i < numAddedObstacles; ++i) {
                int idx = 4 * (i + numObstacles);
                obstacle_renderer.obstacleVa[idx].texCoords = {0.f, 0.f};
                obstacle_renderer.obstacleVa[idx + 1].texCoords = {obstacle_renderer.obstacle_texture_size.x, 0.f};
                obstacle_renderer.obstacleVa[idx + 2].texCoords = {obstacle_renderer.obstacle_texture_size.x, obstacle_renderer.obstacle_texture_size.y};
                obstacle_renderer.obstacleVa[idx + 3].texCoords = {0.f, obstacle_renderer.obstacle_texture_size.y};

                obstacle_renderer.obstacleVa[idx].color = gray;
                obstacle_renderer.obstacleVa[idx + 1].color = gray;
                obstacle_renderer.obstacleVa[idx + 2].color = gray;
                obstacle_renderer.obstacleVa[idx + 3].color = gray;

                auto cellCoords = fluid_attributes.obstaclePositions[numObstacles + i];
                auto pos = fluid_attributes.gridCellToPos(cellCoords.x * fluid_attributes.n + cellCoords.y);
                float px = pos.x;
                float py = pos.y;

                obstacle_renderer.obstacleVa[idx].position = {px - fluid_attributes.halfSpacing, py - fluid_attributes.halfSpacing};
                obstacle_renderer.obstacleVa[idx + 1].position = {px + fluid_attributes.halfSpacing, py - fluid_attributes.halfSpacing};
                obstacle_renderer.obstacleVa[idx + 2].position = {px + fluid_attributes.halfSpacing, py + fluid_attributes.halfSpacing};
                obstacle_renderer.obstacleVa[idx + 3].position = {px - fluid_attributes.halfSpacing, py + fluid_attributes.halfSpacing};
            }
        }  
    }

    void eraseSolids() {
        int localX = static_cast<int>(fluid_attributes.frame_context.world_mouse_pos.x / fluid_attributes.cellSpacing);
        int localY = static_cast<int>(fluid_attributes.frame_context.world_mouse_pos.y / fluid_attributes.cellSpacing);

        int numObstacles = fluid_attributes.obstaclePositions.size();
        int numFreedCells = 0;
        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                int x = localX + i;
                int y = localY + j;

                int idx = x * fluid_attributes.n + y;
                if ((x >= 0 && y >= 0 && x <= fluid_attributes.numX - 1 && y <= fluid_attributes.numY - 1) && (fluid_attributes.cellType[idx] == fluid_attributes.SOLID_CELL && x > 0 && y > 0 && x < fluid_attributes.numX - 1 && y < fluid_attributes.numY - 1)) {
                    fluid_attributes.cellType[idx] = fluid_attributes.AIR_CELL;
                    // find the index of obstacleSet with the position we want to remove
                    // replace that index with the [numObstacles - (numFreedCells + 1)] index of obstacleSet
                    // then just resize the obstacleSet and obstacleVa after this double loop
                    numFreedCells++;

                    int obstacleIdx = find(&fluid_attributes.obstaclePositions, sf::Vector2i{x, y});
                    int vaIdx = 4 * obstacleIdx;

                    int posMv = numObstacles - numFreedCells;
                    int vaMv = 4 * posMv;

                    fluid_attributes.obstaclePositions[obstacleIdx] = fluid_attributes.obstaclePositions[posMv];

                    obstacle_renderer.obstacleVa[vaIdx].position = obstacle_renderer.obstacleVa[vaMv].position;
                    obstacle_renderer.obstacleVa[vaIdx + 1].position = obstacle_renderer.obstacleVa[vaMv + 1].position;
                    obstacle_renderer.obstacleVa[vaIdx + 2].position = obstacle_renderer.obstacleVa[vaMv + 2].position;
                    obstacle_renderer.obstacleVa[vaIdx + 3].position = obstacle_renderer.obstacleVa[vaMv + 3].position;
                }
            }
        }
        obstacle_renderer.obstacleVa.resize(4 * (numObstacles - numFreedCells));
        fluid_attributes.obstaclePositions.resize(numObstacles - numFreedCells);
    }

    void render_objects() {
        if (getPencilActive()) {
            obstacle_renderer.drawPencil(pencilRadius);
        }
    }

    bool getPencilActive() {
        return this->solidDrawing;
    }

    int getPencilRadius() {
        return this->pencilRadius;
    }

    void addToPencilRadius(int add) {
        this->pencilRadius += add;
    }

    void setSolidDrawer(bool set) {
        this->solidDrawing = set;
    }

};