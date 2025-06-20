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

    sf::Clock deltaClock;
    sf::Font font;
    sf::Text text;
    std::ostringstream flipRatio_OSS;
    std::ostringstream gravityX_OSS;
    std::ostringstream gravityY_OSS;
    std::ostringstream k_OSS;
    std::ostringstream vorticity_OSS;
    std::ostringstream pressuerIters_OSS;
    std::ostringstream numParticles_OSS;

    int frame = 0;
    float trueDT;
    float setDT = 1.f / 120.f;
    int fps;

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

        
        // initialize window
        // ----------------------------------------------------------------------------------------------------------------------------
        font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");

        
        text.setFont(font);
        text.setPosition(10, 10);
        text.setFillColor(sf::Color::White);

        
        flipRatio_OSS << std::fixed << std::setprecision(2) << fluid_attributes.flipRatio; 

        gravityX_OSS << std::fixed << std::setprecision(0) << fluid_attributes.gravityX; 

        k_OSS << std::fixed << std::setprecision(1) << 10; 

        vorticity_OSS << std::fixed << std::setprecision(1) << fluid_attributes.vorticityStrength; 

        pressuerIters_OSS << std::fixed << std::setprecision(0) << 30; 

        gravityY_OSS << std::fixed << std::setprecision(0) << fluid_attributes.gravityY; 

        numParticles_OSS << std::fixed << std::setprecision(0) << fluid_attributes.num_particles; 

        window.setFramerateLimit(120);
        // ----------------------------------------------------------------------------------------------------------------------------
    }

    void simulate() {
        window.clear();
        
        sf::Time deltaTime = deltaClock.restart();
        trueDT = deltaTime.asSeconds();

        sf::Vector2i mouse_position = sf::Mouse::getPosition(window);

        fluid_attributes.frame_context.screen_mouse_pos = sf::Vector2f  {static_cast<float>(mouse_position.x), static_cast<float>(mouse_position.y)};
        fluid_attributes.frame_context.world_mouse_pos = scene_renderer.screenToWorld(fluid_attributes.frame_context.screen_mouse_pos);

        fluid_attributes.frame_context.simulation_mouse_pos = fluid_attributes. frame_context.screen_mouse_pos;
        fluid_attributes.frame_context.simulation_mouse_pos = scene_renderer.screenToWorld(fluid_attributes.frame_context.screen_mouse_pos);

        fluid_attributes.frame_context.dt = setDT;

        track_events();
            
        if (!fluid_attributes.stop || fluid_attributes.step) {
            update_environment();
            fluid_attributes.step = false;
        }

        scene_renderer.render_scene();

        handle_zoom();

        render_objects();

        display();
    }

    void update_environment() {
        ++steps;

        // order of need of implementation/optimization:
            // 1) move event handling into scene_handler
            // 2) add a way to switch between particle view and liquid glass view
            // 3) Thread buffers EVERYWHERE
            // 4) implement Multi Grid, then MGPCG
            // 5) make a sampleVelocity(point) function so that you can do RK2 advection easier
            // 6) DDA raycasting & making sure that particle-particle collisions dont push particles into obstacles
            // 7) implement implicit density projection
            // 8) finish cleaning up all this code
            // 9) make it so that you pass in static arrays instead of just numbers of particles in the main file
        
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
            fluid_handler.includeDragObject();
        }

        if (fluid_attributes.vorticityStrength != 0) {
            fluid_handler.applyVorticityConfinementRedBlack();
        }
        
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(miscellaneousTime, duration.count(), steps);

        


        //start = std::chrono::high_resolution_clock::now();
        
        fluid_handler.pressure_solver.projectRedBlackSORMulti(fluid_handler.pressure_solver.numPressureIters);

        // PCG with RBSOR as a preconditioner makes super cool high frequency patterns
        /*fluid_handler.pressure_solver.projectPCG([&]() {
            fluid_handler.pressure_solver.projectRedBlackSORMulti(3); // around 2-6 iterations for preconditioner
        });*/
        //fluid_handler.pressure_solver.projectCG();

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ProjectionTime, duration.count(), steps);



        start = std::chrono::high_resolution_clock::now();
        
        fluid_handler.transfer_grid.TransferToParticles();
        
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        addValueToAverage(ToParticlesTime, duration.count(), steps);
    }

    void track_events() {
        //auto start = std::chrono::high_resolution_clock::now();

        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed) {
            window.close();
            }
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::B) {
                    if (fluid_attributes.getFlipRatio() < 0.99f) {
                        fluid_attributes.addToFlipRatio(0.01);
                        flipRatio_OSS.str("");  
                        flipRatio_OSS.clear();
                        flipRatio_OSS << std::fixed << std::setprecision(2) << fluid_attributes.getFlipRatio(); 
                    }
                }
                else if (event.key.code == sf::Keyboard::S) {
                    if (fluid_attributes.getFlipRatio() > 0.01) {
                        fluid_attributes.addToFlipRatio(-0.01);
                        flipRatio_OSS.str("");  
                        flipRatio_OSS.clear();
                        flipRatio_OSS << std::fixed << std::setprecision(2) << fluid_attributes.getFlipRatio();
                    }
                }
                else if (event.key.code == sf::Keyboard::E) {
                    fluid_attributes.addToVorticityStrength(10);
                    vorticity_OSS.str("");  
                    vorticity_OSS.clear();
                    vorticity_OSS << std::fixed << std::setprecision(1) << fluid_attributes.getVorticityStrength(); 
                }
                else if (event.key.code == sf::Keyboard::W) {
                    if (fluid_attributes.getVorticityStrength() - 10 >= 0.f) {
                        fluid_attributes.addToVorticityStrength(-10);
                        vorticity_OSS.str("");  
                        vorticity_OSS.clear();
                        vorticity_OSS << std::fixed << std::setprecision(1) << fluid_attributes.getVorticityStrength();
                    }
                }
                else if (event.key.code == sf::Keyboard::P) {
                    fluid_handler.pressure_solver.addToNumPressureIters(1);
                    pressuerIters_OSS.str("");  
                    pressuerIters_OSS.clear();
                    pressuerIters_OSS << std::fixed << std::setprecision(0) << fluid_handler.pressure_solver.getNumPressureIters(); 
                }
                else if (event.key.code == sf::Keyboard::O) {
                    if (fluid_handler.pressure_solver.getNumPressureIters() > 0) {
                        fluid_handler.pressure_solver.addToNumPressureIters(-1);
                        pressuerIters_OSS.str("");  
                        pressuerIters_OSS.clear();
                        pressuerIters_OSS << std::fixed << std::setprecision(0) << fluid_handler.pressure_solver.getNumPressureIters();
                    }
                }
                else if (event.key.code == sf::Keyboard::R) {
                    scene_renderer.reset_zoom();
                }
                else if (event.key.code == sf::Keyboard::Num1) {
                    fluid_handler.setDragObjectActive(true);
                    fluid_handler.setForceObjectActive(false);
                    fluid_handler.setGeneratorActive(false);
                    obstacle_handler.setSolidDrawer(false);
                    scene_renderer.setZoomObjectActive(false);
                }
                else if (event.key.code == sf::Keyboard::Num2) {
                    fluid_handler.setDragObjectActive(false);
                    fluid_handler.setForceObjectActive(true);
                    fluid_handler.setGeneratorActive(false);
                    obstacle_handler.setSolidDrawer(false);
                    scene_renderer.setZoomObjectActive(false);
                }
                else if (event.key.code == sf::Keyboard::Num3) {
                    fluid_handler.setDragObjectActive(false);
                    fluid_handler.setForceObjectActive(false);
                    fluid_handler.setGeneratorActive(true);
                    obstacle_handler.setSolidDrawer(false);
                    scene_renderer.setZoomObjectActive(false);
                }
                else if (event.key.code == sf::Keyboard::Num4) {
                    fluid_handler.setDragObjectActive(false);
                    fluid_handler.setForceObjectActive(false);
                    fluid_handler.setGeneratorActive(false);
                    obstacle_handler.setSolidDrawer(true);
                    scene_renderer.setZoomObjectActive(false);
                }
                else if (event.key.code == sf::Keyboard::Num5) {
                    fluid_handler.setDragObjectActive(false);
                    fluid_handler.setForceObjectActive(false);
                    fluid_handler.setGeneratorActive(false);
                    obstacle_handler.setSolidDrawer(false);
                    scene_renderer.setZoomObjectActive(true);
                }
                else if (event.key.code == sf::Keyboard::G) {
                    fluid_attributes.addToGravityY(100);
                    gravityY_OSS.str("");  
                    gravityY_OSS.clear();
                    gravityY_OSS << std::fixed << std::setprecision(0) << fluid_attributes.getGravityY();
                }
                else if (event.key.code == sf::Keyboard::N) {
                    fluid_attributes.addToGravityY(-100);
                    gravityY_OSS.str("");  
                    gravityY_OSS.clear();
                    gravityY_OSS << std::fixed << std::setprecision(0) << fluid_attributes.getGravityY();
                }
                else if (event.key.code == sf::Keyboard::M) {
                    fluid_attributes.addToGravityX(100);
                    gravityX_OSS.str("");  
                    gravityX_OSS.clear();
                    gravityX_OSS << std::fixed << std::setprecision(0) << fluid_attributes.getGravityX();
                }
                else if (event.key.code == sf::Keyboard::H) {
                    fluid_attributes.addToGravityX(-100);
                    gravityX_OSS.str("");  
                    gravityX_OSS.clear();
                    gravityX_OSS << std::fixed << std::setprecision(0) << fluid_attributes.getGravityX();
                }
                else if (event.key.code == sf::Keyboard::C) {
                    if (fluid_handler.pressure_solver.getDivergenceModifier() > 0) {
                        fluid_handler.pressure_solver.addToDivergenceModifier(-1);
                    }
                    k_OSS.str("");  
                    k_OSS.clear();
                    k_OSS << std::fixed << std::setprecision(1) << fluid_handler.pressure_solver.getDivergenceModifier();
                }
                else if (event.key.code == sf::Keyboard::D) {
                    fluid_handler.pressure_solver.addToDivergenceModifier(1);
                    k_OSS.str("");  
                    k_OSS.clear();
                    k_OSS << std::fixed << std::setprecision(1) << fluid_handler.pressure_solver.getDivergenceModifier();
                }
                else if (event.key.code == sf::Keyboard::A) {
                    scene_renderer.fluid_renderer.setNextRenderPattern();
                }
                else if (event.key.code == sf::Keyboard::F) {
                    fluid_attributes.setFireActive(!fluid_attributes.getFireActive());
                }
                else if (event.key.code == sf::Keyboard::Q) {
                    std::cout << 
                    "Fill Grid: " << getFillGridTime() << "\n" <<
                    "Miscellaneous: " << getMiscellaneousTime() << "\n" <<
                    "Collision: " << getCollisionTime() << "\n" <<
                    "Obstacle Collision: " << getObstacleCollisionTime() << "\n" <<
                    "To Grid: " << getToGridTime() << "\n" <<
                    "Density Update: " << getDensityUpdateTime() << "\n" <<
                    "Projection: " << getProjectionTime() << "\n" <<
                    "To Particles: " << getToParticlesTime() << "\n" <<
                    "Rendering: " << getRenderingTime() << "\n";/* <<
                    "Whole Step: " << fluid_handler.getSimStepTime() << "\n" <<
                    "Combined: " << fluid_handler.getCombinedTime() << "\n" <<
                    "Before Sim Step: " << beforeSimStep << "\n" <<
                    "After Sim Step: " << afterSimStep << "\n";*/
                    window.close();
                }
                else if (event.key.code == sf::Keyboard::Y) {
                    fluid_attributes.setStop(!fluid_attributes.getStop());
                }
                else if (event.key.code == sf::Keyboard::U) {
                    fluid_attributes.setStep(true);
                }
            }
            else if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    fluid_attributes.frame_context.leftMouseDown = true;
                    fluid_attributes.frame_context.justPressed = true;
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    fluid_attributes.frame_context.rightMouseDown = true;
                }
            }
            else if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    fluid_attributes.frame_context.leftMouseDown = false;
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    fluid_attributes.frame_context.rightMouseDown = false;
                }
            }
            else if (event.type == sf::Event::MouseWheelScrolled) {

                float mouseWheelDelta = event.mouseWheelScroll.delta;

                if (fluid_handler.getDragObjectActive()) {
                    fluid_handler.addToDragObjectRadius(20 * mouseWheelDelta);
                }

                else if (fluid_handler.getForceObjectActive()) {
                    fluid_handler.addToForceObjectRadius(20 * mouseWheelDelta);
                }

                else if (fluid_handler.getGeneratorActive()) {
                    fluid_handler.addToGeneratorRadius(20 * mouseWheelDelta);
                }

                else if (obstacle_handler.getPencilActive()) {
                    int pencilRadius = obstacle_handler.getPencilRadius();
                    if (!(pencilRadius > fluid_attributes.getNumX() / 2 && mouseWheelDelta > 0)) {
                        obstacle_handler.addToPencilRadius(mouseWheelDelta);
                        if (obstacle_handler.getPencilRadius() < 0) {
                            obstacle_handler.addToPencilRadius(-obstacle_handler.getPencilRadius());
                        }
                    }
                }

                else if (scene_renderer.getZoomObjectActive()) {
                    scene_renderer.wheelZoom(mouseWheelDelta);
                }
            }
            /*end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            addValueToAverage(afterSimStep, duration.count(), numDT);*/
       }
    }

    void display() {
        if (!fluid_attributes.getStop()) {
            frame++;
            if (frame == 20) {
                fps = 1 / trueDT;
                frame = 0;
            }
        }
        else {
            text.setFillColor(sf::Color::Red);
        }

        text.setPosition(fluid_attributes.frame_context.WIDTH - 70, 10);
        text.setString(std::to_string(fps)); 
        window.draw(text);


        if (fluid_attributes.getStop()) {
            text.setFillColor(sf::Color::White);
        }

        text.setPosition(fluid_attributes.frame_context.WIDTH - 175, 10);
        text.setString(flipRatio_OSS.str());
        window.draw(text);

        text.setPosition(fluid_attributes.frame_context.WIDTH - 275, 10);
        text.setString(gravityY_OSS.str());
        window.draw(text);

        text.setPosition(fluid_attributes.frame_context.WIDTH - 375, 10);
        text.setString(gravityX_OSS.str());
        window.draw(text);

        text.setPosition(fluid_attributes.frame_context.WIDTH - 475, 10);
        text.setString(k_OSS.str());
        window.draw(text);

        text.setPosition(fluid_attributes.frame_context.WIDTH - 575, 10);
        text.setString(pressuerIters_OSS.str());
        window.draw(text);

        text.setPosition(fluid_attributes.frame_context.WIDTH - 675, 10);
        text.setString(vorticity_OSS.str());
        window.draw(text);

        if (fluid_handler.getGeneratorActive() && (fluid_attributes.frame_context.leftMouseDown || fluid_attributes.frame_context.rightMouseDown)) {
            numParticles_OSS.str("");  
            numParticles_OSS.clear();
            numParticles_OSS << std::fixed << std::setprecision(1) << fluid_attributes.getNumParticles();
        }

        text.setPosition(fluid_attributes.frame_context.WIDTH - 775, 10);
        text.setString(numParticles_OSS.str());
        window.draw(text);

        window.display();

        fluid_attributes.frame_context.zooming = false;
        fluid_attributes.frame_context.justPressed = false;
    }

    void handle_zoom() {
        if (fluid_attributes.frame_context.zooming) {
            float zoom = fluid_attributes.frame_context.zoom_amount;
            float zoomed_particle_rad = fluid_attributes.radius / zoom;
            fluid_handler.dragObjectSimRadius = fluid_handler.dragObjectRenderRadius / zoom;
            fluid_handler.forceObjectSimRadius = fluid_handler.forceObjectRenderRadius / zoom;
            fluid_handler.checkForceObjectseparationDist = (zoomed_particle_rad + fluid_handler.forceObjectSimRadius) * (zoomed_particle_rad + fluid_handler.forceObjectSimRadius);
            fluid_handler.generatorObjectSimRadius = fluid_handler.generatorObjectRenderRadius / zoom;
        }
    }

    void render_objects() {
        fluid_handler.render_objects(window);
        obstacle_handler.render_objects();
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