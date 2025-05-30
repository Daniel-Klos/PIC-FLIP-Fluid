#include <SFML/Graphics.hpp>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "fluid_handler.hpp"

void addValueToAverage(float& value, float newValue, int steps) {
    value += (newValue - value) / steps;
}

int main()
{
    // Adjust the sim to a good size for you. Note that if you change these numbers, you will have to change the settings
    int WIDTH = 2500; // 2500 
    int HEIGHT = 1300; // 1300

    // SETTINGS
    // ---------------------------------------
    // numparticles and gravity are self explanatory
    // divergenceModifier is a constant which improves convergence of SOR (pressure projection) grid solve; dont put too high or low
    // numpressureiters is how many times SOR iterates
    // diffusionratio is for rendering with the blue and white particles
    // seperationinit is the starting seperation for the particles. make sure that no particles are intersecting and no particles are outside of the bounds of the sim upon initialization. you want them to be evenly spaced out at the start so that a good density sample can be taken
    // vorticitystrength is how strong vorticity confinement forces are if you choose to include that in the sim 

    // multithread
    int numParticles = 5000; 
    float gravityX = 5500.f; // 5500
    float gravityY = 0.f;
    float divergenceModifier = 10.f; // 8 
    int gridNumX = 225; // 100, 150, 200
    int numPressureIters = 30; // 40 
    float diffusionRatio = 0.75f; 
    float flipRatio = 0.9f;
    float vorticityStrength = 250.f;

    // single thread
    /*int numParticles = 5000; 
    float gravityX = 5500.f; // 5500
    float gravityY = 0.f;
    float divergenceModifier = 10.f; // 8 
    int gridNumX = 125; // 100, 150, 200
    int numPressureIters = 30; // 40 
    float diffusionRatio = 0.75f; 
    float flipRatio = 0.9f;
    float vorticityStrength = 250.f;*/

    sf::Font font;
    font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");

    sf::Text text;
    text.setFont(font);
    text.setPosition(10, 10);
    text.setFillColor(sf::Color::White);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << flipRatio; 

    std::ostringstream oss2;
    oss2 << std::fixed << std::setprecision(0) << gravityX; 

    std::ostringstream oss3;
    oss3 << std::fixed << std::setprecision(1) << divergenceModifier; 

    std::ostringstream oss4;
    oss4 << std::fixed << std::setprecision(1) << vorticityStrength; 

    std::ostringstream oss5;
    oss5 << std::fixed << std::setprecision(0) << numPressureIters; 

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

    const uint32_t numThreads = std::min(static_cast<uint32_t>(16), maxThreads); // 11

    tp::ThreadPool thread_pool(numThreads);

    const float overRelaxation = 1.9f; // 1.9

    FluidState fluid_attributes = FluidState(numParticles, WIDTH, HEIGHT, gridNumX, vorticityStrength, flipRatio, gravityX, gravityY, thread_pool);

    PressureSolver pressure_solver = PressureSolver(fluid_attributes, numPressureIters);

    TransferGrid transfer_grid = TransferGrid(fluid_attributes);

    FluidRenderer fluid_renderer = FluidRenderer(fluid_attributes, window);

    FluidHandler fluid = FluidHandler(divergenceModifier, overRelaxation, numPressureIters, fluid_attributes, pressure_solver, transfer_grid, fluid_renderer);

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
                        oss.str("");  
                        oss.clear();
                        oss << std::fixed << std::setprecision(2) << fluid_attributes.getFlipRatio(); 
                    }
                }
                else if (event.key.code == sf::Keyboard::S) {
                    if (fluid_attributes.getFlipRatio() > 0.01) {
                        fluid_attributes.addToFlipRatio(-0.01);
                        oss.str("");  
                        oss.clear();
                        oss << std::fixed << std::setprecision(2) << fluid_attributes.getFlipRatio();
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
                    pressure_solver.addToNumPressureIters(1);
                    oss5.str("");  
                    oss5.clear();
                    oss5 << std::fixed << std::setprecision(0) << pressure_solver.getNumPressureIters(); 
                }
                else if (event.key.code == sf::Keyboard::O) {
                    if (pressure_solver.getNumPressureIters() > 0) {
                        pressure_solver.addToNumPressureIters(-1);
                        oss5.str("");  
                        oss5.clear();
                        oss5 << std::fixed << std::setprecision(0) << pressure_solver.getNumPressureIters();
                    }
                }
                else if (event.key.code == sf::Keyboard::Num1) {
                    fluid.setRigidObjectActive(true);
                    fluid.setForceObjectActive(false);
                    fluid.setGeneratorActive(false);
                    fluid.setSolidDrawer(false);
                }
                else if (event.key.code == sf::Keyboard::Num2) {
                    fluid.setRigidObjectActive(false);
                    fluid.setForceObjectActive(true);
                    fluid.setGeneratorActive(false);
                    fluid.setSolidDrawer(false);
                }
                else if (event.key.code == sf::Keyboard::Num3) {
                    fluid.setRigidObjectActive(false);
                    fluid.setForceObjectActive(false);
                    fluid.setGeneratorActive(true);
                    fluid.setSolidDrawer(false);
                }
                else if (event.key.code == sf::Keyboard::Num4) {
                    fluid.setRigidObjectActive(false);
                    fluid.setForceObjectActive(false);
                    fluid.setGeneratorActive(false);
                    fluid.setSolidDrawer(true);
                }
                else if (event.key.code == sf::Keyboard::G) {
                    fluid_attributes.addToGravityX(100);
                    oss2.str("");  
                    oss2.clear();
                    oss2 << std::fixed << std::setprecision(0) << fluid_attributes.getGravityX();
                }
                else if (event.key.code == sf::Keyboard::N) {
                    fluid_attributes.addToGravityX(-100);
                    oss2.str("");  
                    oss2.clear();
                    oss2 << std::fixed << std::setprecision(0) << fluid_attributes.getGravityX();
                }
                else if (event.key.code == sf::Keyboard::M) {
                    fluid_attributes.addToGravityY(100);
                    oss6.str("");  
                    oss6.clear();
                    oss6 << std::fixed << std::setprecision(0) << fluid_attributes.getGravityY();
                }
                else if (event.key.code == sf::Keyboard::H) {
                    fluid_attributes.addToGravityY(-100);
                    oss6.str("");  
                    oss6.clear();
                    oss6 << std::fixed << std::setprecision(0) << fluid_attributes.getGravityY();
                }
                else if (event.key.code == sf::Keyboard::C) {
                    if (pressure_solver.getDivergenceModifier() > 0) {
                        pressure_solver.addToDivergenceModifier(-1);
                    }
                    oss3.str("");  
                    oss3.clear();
                    oss3 << std::fixed << std::setprecision(1) << pressure_solver.getDivergenceModifier();
                }
                else if (event.key.code == sf::Keyboard::D) {
                    pressure_solver.addToDivergenceModifier(1);
                    oss3.str("");  
                    oss3.clear();
                    oss3 << std::fixed << std::setprecision(1) << pressure_solver.getDivergenceModifier();
                }
                else if (event.key.code == sf::Keyboard::A) {
                    fluid.setNextRenderPattern();
                }
                else if (event.key.code == sf::Keyboard::F) {
                    fluid_attributes.setFireActive(!fluid_attributes.getFireActive());
                }
                else if (event.key.code == sf::Keyboard::Q) {
                    /*std::cout << 
                    /*"Fill Grid: " << fluid.getFillGridTime() << "\n" <<
                    "Miscellaneous: " << fluid.getMiscellaneousTime() << "\n" <<
                    "Collision: " << fluid.getCollisionTime() << "\n" <<
                    "Obstacle Collision: " << fluid.getObstacleCollisionTime() << "\n" <<*/
                    /*"To Grid: " << fluid.getToGridTime() << "\n";// <<*/
                    /*"Density Update: " << fluid.getDensityUpdateTime() << "\n" <<
                    "Projection: " << fluid.getProjectionTime() << "\n" <<
                    "To Particles: " << fluid.getToParticlesTime() << "\n" <<
                    "Rendering: " << fluid.getRenderingTime() << "\n" <<*/
                    /*"Whole Step: " << fluid.getSimStepTime() << "\n"; <<
                    "Combined: " << fluid.getCombinedTime() << "\n" <<
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

                if (fluid.getForceObjectActive()) {
                    fluid.addToForceObjectRadius(20 * mouseWheelDelta);
                }

                else if (fluid.getGeneratorActive()) {
                    fluid.addToGeneratorRadius(20 * mouseWheelDelta);
                }

                else if (fluid.getPencilActive()) {
                    int pencilRadius = fluid.getPencilRadius();
                    if (!(pencilRadius > fluid.getNumX() / 2 && mouseWheelDelta > 0)) {
                        fluid.addToPencilRadius(mouseWheelDelta);
                        if (fluid.getPencilRadius() < 0) {
                            fluid.addToPencilRadius(-fluid.getPencilRadius());
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
        fluid.update(setDT, window, leftMouseDown, rightMouseDown, justPressed);
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
        text.setString(oss.str());
        window.draw(text);

        text.setPosition(WIDTH - 275, 10);
        text.setString(oss2.str());
        window.draw(text);

        text.setPosition(WIDTH - 350, 10);
        text.setString(oss6.str());
        window.draw(text);

        text.setPosition(WIDTH - 475, 10);
        text.setString(oss3.str());
        window.draw(text);

        text.setPosition(WIDTH - 575, 10);
        text.setString(oss4.str());
        window.draw(text);

        text.setPosition(WIDTH - 675, 10);
        text.setString(oss5.str());
        window.draw(text);

        if (fluid.getGeneratorActive() && (leftMouseDown || rightMouseDown)) {
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
