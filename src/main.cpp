#include <SFML/Graphics.hpp>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "scene_handler.hpp"

int main()
{
    // Adjust the sim to a good size for you. Note that if you change these numbers, you will have to change the settings
    int WIDTH = 2500; // 2500 
    int HEIGHT = 1300; // 1300

    // SETTINGS
    // ---------------------------------------
    // numparticles and gravity are self explanatory
    // diffusionratio is for rendering with the blue and white particles
    // seperationinit is the starting seperation for the particles. make sure that no particles are intersecting and no particles are outside of the bounds of the sim upon initialization. you want them to be evenly spaced out at the start so that a good density sample can be taken
    // vorticitystrength is how strong vorticity confinement forces are if you choose to include that in the sim 

    // multithread -- low power mode
    int numParticles = 20000; // 20000 --- start with a kinda large number so that good density sample is taken at the start of the simulation
    float gravityY = 5500.f; // 5500
    float gravityX = 0.f;
    int gridNumX = 350; // 350
    float diffusionRatio = 0.75f; 
    float flipRatio = 0.9f;
    float vorticityStrength = 0.f;
    uint32_t numThreads = 6; // 6

    // single thread
    /*int numParticles = 5000; 
    float gravityY = 5500.f; // 5500
    float gravityX = 0.f;
    float divergenceModifier = 10.f; // 8 
    int gridNumX = 300; 
    int numPressureIters = 30; // 40 
    float diffusionRatio = 0.75f; 
    float flipRatio = 0.9f;
    float vorticityStrength = 250.f;
    uint32_t numThreads = 1; // 6*/

    sf::Font font;
    font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");

    sf::Text text;
    text.setFont(font);
    text.setPosition(10, 10);
    text.setFillColor(sf::Color::White);

    std::ostringstream oss1;
    oss1 << std::fixed << std::setprecision(2) << flipRatio; 

    std::ostringstream oss2;
    oss2 << std::fixed << std::setprecision(0) << gravityX; 

    std::ostringstream oss3;
    oss3 << std::fixed << std::setprecision(1) << 10; 

    std::ostringstream oss4;
    oss4 << std::fixed << std::setprecision(1) << vorticityStrength; 

    std::ostringstream oss5;
    oss5 << std::fixed << std::setprecision(0) << 30; 

    std::ostringstream oss6;
    oss6 << std::fixed << std::setprecision(0) << gravityY; 

    std::ostringstream oss7;
    oss7 << std::fixed << std::setprecision(0) << numParticles; 

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "FLIP Simulation");

    sf::Clock deltaClock;

    window.setFramerateLimit(120); // 120

    int frame = 0;
    int fps = 0;

    float interactionRadius = 100.f;

    float interactionStrength;

    bool leftMouseDown = false;
    bool rightMouseDown = false;

    const uint32_t maxThreads = std::thread::hardware_concurrency();

    numThreads = std::min(numThreads, maxThreads); // 10, 11, 16

    tp::ThreadPool thread_pool(numThreads);

    const float overRelaxation = 1.9f; // 1.9

    FrameContext frame_context = FrameContext(WIDTH, HEIGHT);

    FluidState fluid_attributes = FluidState(numParticles, gridNumX, vorticityStrength, flipRatio, gravityX, gravityY, frame_context, thread_pool);

    SceneHandler scene_handler = SceneHandler(fluid_attributes, window);

    bool justPressed = false;

    bool forceObjectActive = true; 

    float totalDT = 0;
    float numDT = 0;

    float afterSimStep = 0.f;
    float beforeSimStep = 0.f;

    while (window.isOpen())
    {

        //auto start = std::chrono::high_resolution_clock::now();

        sf::Time deltaTime = deltaClock.restart();
        float setDT = 1.f / 120.f;
        float trueDT = deltaTime.asSeconds();

        //totalDT += dt;
        numDT++;

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
                        oss1.str("");  
                        oss1.clear();
                        oss1 << std::fixed << std::setprecision(2) << fluid_attributes.getFlipRatio(); 
                    }
                }
                else if (event.key.code == sf::Keyboard::S) {
                    if (fluid_attributes.getFlipRatio() > 0.01) {
                        fluid_attributes.addToFlipRatio(-0.01);
                        oss1.str("");  
                        oss1.clear();
                        oss1 << std::fixed << std::setprecision(2) << fluid_attributes.getFlipRatio();
                    }
                }
                else if (event.key.code == sf::Keyboard::E) {
                    fluid_attributes.addToVorticityStrength(10);
                    oss4.str("");  
                    oss4.clear();
                    oss4 << std::fixed << std::setprecision(1) << fluid_attributes.getVorticityStrength(); 
                }
                else if (event.key.code == sf::Keyboard::W) {
                    if (fluid_attributes.getVorticityStrength() - 10 >= 0.f) {
                        fluid_attributes.addToVorticityStrength(-10);
                        oss4.str("");  
                        oss4.clear();
                        oss4 << std::fixed << std::setprecision(1) << fluid_attributes.getVorticityStrength();
                    }
                }
                else if (event.key.code == sf::Keyboard::P) {
                    scene_handler.fluid_handler.pressure_solver.addToNumPressureIters(1);
                    oss5.str("");  
                    oss5.clear();
                    oss5 << std::fixed << std::setprecision(0) << scene_handler.fluid_handler.pressure_solver.getNumPressureIters(); 
                }
                else if (event.key.code == sf::Keyboard::O) {
                    if (scene_handler.fluid_handler.pressure_solver.getNumPressureIters() > 0) {
                        scene_handler.fluid_handler.pressure_solver.addToNumPressureIters(-1);
                        oss5.str("");  
                        oss5.clear();
                        oss5 << std::fixed << std::setprecision(0) << scene_handler.fluid_handler.pressure_solver.getNumPressureIters();
                    }
                }
                else if (event.key.code == sf::Keyboard::Num1) {
                    scene_handler.fluid_handler.setRigidObjectActive(true);
                    scene_handler.fluid_handler.setForceObjectActive(false);
                    scene_handler.fluid_handler.setGeneratorActive(false);
                    scene_handler.obstacle_handler.setSolidDrawer(false);
                }
                else if (event.key.code == sf::Keyboard::Num2) {
                    scene_handler.fluid_handler.setRigidObjectActive(false);
                    scene_handler.fluid_handler.setForceObjectActive(true);
                    scene_handler.fluid_handler.setGeneratorActive(false);
                    scene_handler.obstacle_handler.setSolidDrawer(false);
                }
                else if (event.key.code == sf::Keyboard::Num3) {
                    scene_handler.fluid_handler.setRigidObjectActive(false);
                    scene_handler.fluid_handler.setForceObjectActive(false);
                    scene_handler.fluid_handler.setGeneratorActive(true);
                    scene_handler.obstacle_handler.setSolidDrawer(false);
                }
                else if (event.key.code == sf::Keyboard::Num4) {
                    scene_handler.fluid_handler.setRigidObjectActive(false);
                    scene_handler.fluid_handler.setForceObjectActive(false);
                    scene_handler.fluid_handler.setGeneratorActive(false);
                    scene_handler.obstacle_handler.setSolidDrawer(true);
                }
                else if (event.key.code == sf::Keyboard::G) {
                    fluid_attributes.addToGravityY(100);
                    oss6.str("");  
                    oss6.clear();
                    oss6 << std::fixed << std::setprecision(0) << fluid_attributes.getGravityY();
                }
                else if (event.key.code == sf::Keyboard::N) {
                    fluid_attributes.addToGravityY(-100);
                    oss6.str("");  
                    oss6.clear();
                    oss6 << std::fixed << std::setprecision(0) << fluid_attributes.getGravityY();
                }
                else if (event.key.code == sf::Keyboard::M) {
                    fluid_attributes.addToGravityX(100);
                    oss2.str("");  
                    oss2.clear();
                    oss2 << std::fixed << std::setprecision(0) << fluid_attributes.getGravityX();
                }
                else if (event.key.code == sf::Keyboard::H) {
                    fluid_attributes.addToGravityX(-100);
                    oss2.str("");  
                    oss2.clear();
                    oss2 << std::fixed << std::setprecision(0) << fluid_attributes.getGravityX();
                }
                else if (event.key.code == sf::Keyboard::C) {
                    if (scene_handler.fluid_handler.pressure_solver.getDivergenceModifier() > 0) {
                        scene_handler.fluid_handler.pressure_solver.addToDivergenceModifier(-1);
                    }
                    oss3.str("");  
                    oss3.clear();
                    oss3 << std::fixed << std::setprecision(1) << scene_handler.fluid_handler.pressure_solver.getDivergenceModifier();
                }
                else if (event.key.code == sf::Keyboard::D) {
                    scene_handler.fluid_handler.pressure_solver.addToDivergenceModifier(1);
                    oss3.str("");  
                    oss3.clear();
                    oss3 << std::fixed << std::setprecision(1) << scene_handler.fluid_handler.pressure_solver.getDivergenceModifier();
                }
                else if (event.key.code == sf::Keyboard::A) {
                    scene_handler.scene_renderer.fluid_renderer.setNextRenderPattern();
                }
                else if (event.key.code == sf::Keyboard::F) {
                    fluid_attributes.setFireActive(!fluid_attributes.getFireActive());
                }
                else if (event.key.code == sf::Keyboard::Q) {
                    std::cout << 
                    "Fill Grid: " << scene_handler.getFillGridTime() << "\n" <<
                    "Miscellaneous: " << scene_handler.getMiscellaneousTime() << "\n" <<
                    "Collision: " << scene_handler.getCollisionTime() << "\n" <<
                    "Obstacle Collision: " << scene_handler.getObstacleCollisionTime() << "\n" <<
                    "To Grid: " << scene_handler.getToGridTime() << "\n" <<
                    "Density Update: " << scene_handler.getDensityUpdateTime() << "\n" <<
                    "Projection: " << scene_handler.getProjectionTime() << "\n" <<
                    "To Particles: " << scene_handler.getToParticlesTime() << "\n" <<
                    "Rendering: " << scene_handler.getRenderingTime() << "\n";/* <<
                    "Whole Step: " << scene_handler.fluid_handler.getSimStepTime() << "\n" <<
                    "Combined: " << scene_handler.fluid_handler.getCombinedTime() << "\n" <<
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
                    leftMouseDown = true;
                    justPressed = true;
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    rightMouseDown = true;
                }
            }
            else if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    leftMouseDown = false;
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    rightMouseDown = false;
                }
            }
            else if (event.type == sf::Event::MouseWheelMoved) {

                int mouseWheelDelta = event.mouseWheel.delta;

                if (scene_handler.fluid_handler.getForceObjectActive()) {
                    scene_handler.fluid_handler.addToForceObjectRadius(20 * mouseWheelDelta);
                }

                else if (scene_handler.fluid_handler.getGeneratorActive()) {
                    scene_handler.fluid_handler.addToGeneratorRadius(20 * mouseWheelDelta);
                }

                else if (scene_handler.obstacle_handler.getPencilActive()) {
                    int pencilRadius = scene_handler.obstacle_handler.getPencilRadius();
                    if (!(pencilRadius > scene_handler.fluid_attributes.getNumX() / 2 && mouseWheelDelta > 0)) {
                        scene_handler.obstacle_handler.addToPencilRadius(mouseWheelDelta);
                        if (scene_handler.obstacle_handler.getPencilRadius() < 0) {
                            scene_handler.obstacle_handler.addToPencilRadius(-scene_handler.obstacle_handler.getPencilRadius());
                        }
                    }
                }
            }
        }

        window.clear();

        /*auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        addValueToAverage(beforeSimStep, duration.count(), numDT);*/

        //auto start = std::chrono::high_resolution_clock::now();
        scene_handler.simulate(setDT, leftMouseDown, rightMouseDown, justPressed);
        //auto end = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        //start = std::chrono::high_resolution_clock::now();
        if (!fluid_attributes.getStop()) {
            frame++;
            if (frame == 20) {
                fps = (int)(1.f / trueDT);
                frame = 0;
            }
        }
        else {
            text.setFillColor(sf::Color::Red);
        }

        text.setPosition(WIDTH - 70, 10);
        text.setString(std::to_string(fps)); 
        window.draw(text);


        if (fluid_attributes.getStop()) {
            text.setFillColor(sf::Color::White);
        }

        text.setPosition(WIDTH - 175, 10);
        text.setString(oss1.str());
        window.draw(text);

        text.setPosition(WIDTH - 275, 10);
        text.setString(oss5.str());
        window.draw(text);

        text.setPosition(WIDTH - 375, 10);
        text.setString(oss3.str());
        window.draw(text);

        text.setPosition(WIDTH - 475, 10);
        text.setString(oss4.str());
        window.draw(text);

        text.setPosition(WIDTH - 575, 10);
        text.setString(oss2.str());
        window.draw(text);

        text.setPosition(WIDTH - 675, 10);
        text.setString(oss6.str());
        window.draw(text);

        if (scene_handler.fluid_handler.getGeneratorActive() && (leftMouseDown || rightMouseDown)) {
            oss7.str("");  
            oss7.clear();
            oss7 << std::fixed << std::setprecision(1) << fluid_attributes.getNumParticles();
        }
        
        text.setPosition(WIDTH - 775, 10);
        text.setString(oss7.str());
        window.draw(text);


        window.display();
 
        justPressed = false;
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        addValueToAverage(afterSimStep, duration.count(), numDT);*/
       
    }

    return 0;
}