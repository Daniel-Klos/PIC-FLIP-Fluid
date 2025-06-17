#pragma once
#include "simulation_state.hpp"
#include "fluid_handler.hpp"
#include "obstacle_handler.hpp"
#include "scene_renderer.hpp"

struct SceneHandler {
    FluidState &fluid_attributes;
    sf::RenderWindow &window;

    SceneRenderer scene_renderer;
    ObstacleHandler obstacle_handler;
    FluidHandler fluid_handler;

    float mouseX;
    float mouseY;

    float steps = 0.f;
    float CollisionTime = 0.f;
    float ObstacleCollisionTime = 0.f;
    float ToGridTime = 0.f;
    float DensityUpdateTime = 0.f;
    float ProjectionTime = 0.f;
    float ToParticlesTime = 0.f;
    float RenderingTime = 0.f;
    float FillGridTime = 0.f;
    float SimStepTime = 0.f;
    float miscellaneousTime = 0.f;

    SceneHandler(FluidState &fas, sf::RenderWindow &w)
        : fluid_attributes(fas), 
          window(w),
          scene_renderer(fas, w), 
          obstacle_handler(fas, scene_renderer.obstacle_renderer), 
          fluid_handler(fas, scene_renderer.fluid_renderer)
    {
        // initialize obstacle positions
        // ----------------------------------------------------------------------------------------------------------------------------

        size_t numInitObstacles = 2 * fluid_attributes.numX + 2 * (fluid_attributes.numY - 2);
        fluid_attributes.obstaclePositions.resize(numInitObstacles);
        int idx = 0;
        for (int i = 0; i < fluid_attributes.numX; ++i) {
            fluid_attributes.cellType[i * fluid_attributes.n] = fluid_attributes.SOLID_CELL;
            fluid_attributes.obstaclePositions[idx] = sf::Vector2i{i, 0};
            ++idx;
        }
        for (int i = 0; i < fluid_attributes.numX; ++i) {
            fluid_attributes.cellType[i * fluid_attributes.n + fluid_attributes.numY - 1] = fluid_attributes.SOLID_CELL;
            fluid_attributes.obstaclePositions[idx] = sf::Vector2i{i, fluid_attributes.numY - 1};
            ++idx;
        }
        for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
            fluid_attributes.cellType[j] = fluid_attributes.SOLID_CELL;
            fluid_attributes.obstaclePositions[idx] = sf::Vector2i{0, j};
            ++idx;
        }
        for (int j = 1; j < fluid_attributes.numY - 1; ++j) {
            fluid_attributes.cellType[fluid_attributes.numX * fluid_attributes.numY + j - fluid_attributes.n] = fluid_attributes.SOLID_CELL;
            fluid_attributes.obstaclePositions[idx] = sf::Vector2i{fluid_attributes.numX - 1, j};
            ++idx;
        }
        // ----------------------------------------------------------------------------------------------------------------------------



        // initialize obstacle_renderer
        // ----------------------------------------------------------------------------------------------------------------------------
    
        scene_renderer.obstacle_renderer.obstacleVa.resize(4 * numInitObstacles);

        sf::Color gray = sf::Color(150, 150, 150);

        scene_renderer.obstacle_renderer.obstacleTexture.loadFromFile("gray_square.png");
        scene_renderer.obstacle_renderer.obstacleTexture.generateMipmap();
        scene_renderer.obstacle_renderer.obstacle_texture_size = static_cast<sf::Vector2f>(scene_renderer.obstacle_renderer.obstacleTexture.getSize());
        for (int index = 0; index < numInitObstacles; ++index) {
            int i = 4 * index;
            scene_renderer.obstacle_renderer.obstacleVa[i].texCoords = {0.f, 0.f};
            scene_renderer.obstacle_renderer.obstacleVa[i + 1].texCoords = {scene_renderer.obstacle_renderer.obstacle_texture_size.x, 0.f};
            scene_renderer.obstacle_renderer.obstacleVa[i + 2].texCoords = {scene_renderer.obstacle_renderer.obstacle_texture_size.x, scene_renderer.obstacle_renderer.obstacle_texture_size.y};
            scene_renderer.obstacle_renderer.obstacleVa[i + 3].texCoords = {0.f, scene_renderer.obstacle_renderer.obstacle_texture_size.y};

            scene_renderer.obstacle_renderer.obstacleVa[i].color = gray;
            scene_renderer.obstacle_renderer.obstacleVa[i + 1].color = gray;
            scene_renderer.obstacle_renderer.obstacleVa[i + 2].color = gray;
            scene_renderer.obstacle_renderer.obstacleVa[i + 3].color = gray;
        }
        scene_renderer.obstacle_renderer.obstacleStates.texture = &scene_renderer.obstacle_renderer.obstacleTexture;
            
        for (int i = 0; i < numInitObstacles; ++i) {
            int gx = fluid_attributes.obstaclePositions[i].x;
            int gy = fluid_attributes.obstaclePositions[i].y;
            auto pos = fluid_attributes.gridCellToPos(gx * fluid_attributes.n + gy);
            float px = pos.x;
            float py = pos.y;
            size_t vaIdx = 4 * i;
            scene_renderer.obstacle_renderer.obstacleVa[vaIdx].position = {px - fluid_attributes.halfSpacing, py - fluid_attributes.halfSpacing};
            scene_renderer.obstacle_renderer.obstacleVa[vaIdx + 1].position = {px + fluid_attributes.halfSpacing, py - fluid_attributes.halfSpacing};
            scene_renderer.obstacle_renderer.obstacleVa[vaIdx + 2].position = {px + fluid_attributes.halfSpacing, py + fluid_attributes.halfSpacing};
            scene_renderer.obstacle_renderer.obstacleVa[vaIdx + 3].position = {px - fluid_attributes.halfSpacing, py + fluid_attributes.halfSpacing};
        }
        // ----------------------------------------------------------------------------------------------------------------------------
    }

    void simulate(float dt_, bool leftMouseDown, bool rightMouseDown, bool justPressed) {
        sf::Vector2i mouse_position = sf::Mouse::getPosition(window);

        fluid_attributes.frame_context.screen_mouse_pos = sf::Vector2f{static_cast<float>(mouse_position.x), static_cast<float>(mouse_position.y)};
        fluid_attributes.frame_context.world_mouse_pos = scene_renderer.screenToWorld(fluid_attributes.frame_context.screen_mouse_pos);
        fluid_attributes.frame_context.simulation_mouse_pos = fluid_attributes.frame_context.screen_mouse_pos;
            
        fluid_attributes.frame_context.simulation_mouse_pos = scene_renderer.screenToWorld(fluid_attributes.frame_context.screen_mouse_pos);

        fluid_attributes.frame_context.leftMouseDown = leftMouseDown;
        fluid_attributes.frame_context.rightMouseDown = rightMouseDown;
        fluid_attributes.frame_context.justPressed = justPressed;

        fluid_attributes.frame_context.dt = dt_;
        
        if (!fluid_attributes.stop || fluid_attributes.step) {
            update_environment();
            fluid_attributes.step = false;
        }

        scene_renderer.render_scene();
        render_objects();
    }

    void update_environment() {
        ++steps;

        // order of need of implementation/optimization:
            // 1) implement zooming in and out (make a zoom object)
                // 1) compute num_in_view (thread buffers) and fill a boolean array of the particle being in view or not AND THEN fill a vector in_view, which can be resized and then filled in parallel. Just insert this logic into existing loops
            // 2) move event handling into scene_handler
            // 2) add a way to switch between particle view and liquid glass view
            // 2) Thread buffers EVERYWHERE
            // 3) finish cleaning up all this code
            // 4) incompressibility -- implement MGPCG
            // 5) make a sampleVelocity(point) function so that you can do RK2 advection easier
            // 6) make it so that you pass in static arrays instead of just numbers of particles in the main file
            // 7) level set & fast sweeping for separation from obstacles, or DDA raycasting & making sure that particle-particle collisions dont push particles into obstacles
            // 8) implement implicit density projection
        
        auto start = std::chrono::high_resolution_clock::now();

        fluid_handler.integrateMulti();

        // fluid_handler.pressure_solver.project_density_implicit();

        if (fluid_handler.fluid_renderer.renderPattern == 3 && fluid_attributes.fireActive) {
            fluid_handler.makeFireMulti();
        }

        if (fluid_attributes.fireActive) {
            std::fill(begin(fluid_handler.collisions), end(fluid_handler.collisions), 0);
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count(), steps);


        start = std::chrono::high_resolution_clock::now();
        fluid_handler.addObjectsToGrids();

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(FillGridTime, duration.count(), steps);

        


        start = std::chrono::high_resolution_clock::now();

        if (obstacle_handler.solidDrawing && fluid_attributes.frame_context.leftMouseDown) {
            obstacle_handler.drawSolids();
        } else if (obstacle_handler.solidDrawing && fluid_attributes.frame_context.rightMouseDown) {
            obstacle_handler.eraseSolids();
        }

        if (fluid_handler.generatorActive && fluid_attributes.frame_context.leftMouseDown) {
            fluid_handler.generate();
        } else if (fluid_handler.generatorActive && fluid_attributes.frame_context.rightMouseDown) {
            fluid_handler.remove();
        }
    
        if (fluid_handler.forceObjectActive && fluid_attributes.frame_context.leftMouseDown) {
            fluid_handler.makeForceObjectQueries(-250); // pulling, -250
        } else if (fluid_handler.forceObjectActive && fluid_attributes.frame_context.rightMouseDown) { 
            fluid_handler.makeForceObjectQueries(1000); // pushing, 1000
        }

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count(), steps);


        

        start = std::chrono::high_resolution_clock::now();
        fluid_handler.solveCollisions();
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(CollisionTime, duration.count(), steps);




        start = std::chrono::high_resolution_clock::now();

        if (fluid_handler.fluid_renderer.renderPattern == 3) {
            fluid_handler.heatGroundMulti();
        }

        obstacle_handler.collideSurfacesMulti();

        obstacle_handler.constrainWallsMulti();

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ObstacleCollisionTime, duration.count(), steps);




        start = std::chrono::high_resolution_clock::now();
    
        fluid_handler.transfer_grid.TransferToGrid();
        
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ToGridTime, duration.count(), steps);

        start = std::chrono::high_resolution_clock::now();
        
        fluid_handler.transfer_grid.updateCellDensitiesMulti();

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(DensityUpdateTime, duration.count(), steps);





        start = std::chrono::high_resolution_clock::now();
        
        if (fluid_handler.dragObjectActive) {
            fluid_handler.includeDragObject(fluid_attributes.frame_context.leftMouseDown, fluid_attributes.frame_context.justPressed);
        }

        if (fluid_attributes.vorticityStrength != 0) {
            fluid_handler.applyVorticityConfinementRedBlack();
        }
        
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count(), steps);

        


        //start = std::chrono::high_resolution_clock::now();
        
        //pressure_solver.projectMICCG(pressure_solver.numPressureIters);
        //pressure_solver.projectSOR(pressure_solver.numPressureIters);
        //pressure_solver.projectRedBlackGS(pressure_solver.numPressureIters);
        fluid_handler.pressure_solver.projectRedBlackGSMulti(fluid_handler.pressure_solver.numPressureIters, fluid_attributes.numThreads);
        //pressure_solver.projectCG(pressure_solver.numPressureIters);

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ProjectionTime, duration.count(), steps);



        start = std::chrono::high_resolution_clock::now();
        
        fluid_handler.transfer_grid.TransferToParticles();
        
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ToParticlesTime, duration.count(), steps);
    }

    void render_objects() {
        fluid_handler.render_objects(window);
        obstacle_handler.render_objects();
    }

    void zoom_objects() {

    }

    float getCombinedTime() {
        return FillGridTime + miscellaneousTime + CollisionTime + ObstacleCollisionTime + ToGridTime + DensityUpdateTime + ProjectionTime + ToParticlesTime + RenderingTime;
    }

    float getSimStepTime() {
        return SimStepTime;
    }

    float getFillGridTime() {
        return FillGridTime;
    }

    float getMiscellaneousTime() {
        return miscellaneousTime;
    }

    float getCollisionTime() {
        return CollisionTime;
    }

    float getObstacleCollisionTime() {
        return ObstacleCollisionTime;
    }

    float getToGridTime() {
        return ToGridTime;
    }

    float getDensityUpdateTime() {
        return DensityUpdateTime;
    }

    float getProjectionTime() {
        return ProjectionTime;
    }

    float getToParticlesTime() {
        return ToParticlesTime;
    }

    float getRenderingTime() {
        return RenderingTime;
    }

};