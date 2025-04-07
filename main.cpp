#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <array>
#include <sstream>
#include <iomanip>
#include <thread>

#include "fluid.hpp"

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

    // experiment
    /*int numParticles = 15000; 
    float gravity = 5500.f; 
    float divergenceModifier = 8.5f; 
    float gridSize = 70.f; 
    int numPressureIters = 50; 
    float diffusionRatio = 0.95f; 
    float flipRatio = 0.9f;
    float seperationInit = 2.7f;
    float vorticityStrength = 10.5f;//19.5f; */

    // space
    /*int numParticles = 1000; 
    float gravity = 0.f; 
    float divergenceModifier = 15.f; 
    float gridSize = 90.f; 
    int numPressureIters = 30; 
    float diffusionRatio = 0.8f;
    float flipRatio = 0.8f;
    float seperationInit = 2.7f;
    float vorticityStrength = 4.f;*/

    // lots
    /*int32_t numParticles = 25000; 
    float gravity = 5400.f; 
    float divergenceModifier = 9.5f; 
    float gridSize = 80.f; 
    int32_t numPressureIters = 30; 
    float diffusionRatio = 0.9f; 
    float flipRatio = 0.9f;
    float seperationInit = 2.4f;
    float vorticityStrength = 0.f;*/

    // fire
    /*int32_t numParticles = 14500; //13000
    float gravity = 5500.f; 
    float divergenceModifier = 40.f; 
    float gridSize = 50.f;  
    int32_t numPressureIters = 30; // 30 
    float diffusionRatio = 0.85f; 
    float flipRatio = 0.9f;
    float seperationInit = 4.13f; 
    float vorticityStrength = 5.f;//19.5f; */

    // grid test
    /*int32_t numParticles = 30000; 
    float gravity = 5500.f; 
    float divergenceModifier = 30.f; // 8.5 for 15000 
    float gridSize = 100.f; // 100
    int32_t numPressureIters = 30; // 30 
    float diffusionRatio = 0.85f; 
    float flipRatio = 0.9f;
    float seperationInit = 1.9f; // 1.9 for 100, 2.4 for 80
    float vorticityStrength = 5.f;//19.5f; */

    // a lot
    /*int numParticles = 26400; 
    float gravity = 4500.f; 
    float divergenceModifier = 10.f;
    float gridSize = 85.f; 
    int numPressureIters = 30;
    float diffusionRatio = 0.9f; 
    float flipRatio = 0.80f;
    float seperationInit = 2.f;
    float vorticityStrength = 150.f;*/

    // put laptop on performance mode
    /*int numParticles = 35000; 
    float gravity = 4500.f; // 4500
    float divergenceModifier = 17.f;
    float gridSize = 90.f; 
    int numPressureIters = 30;
    float diffusionRatio = 0.9f; 
    float flipRatio = 0.80f;
    float seperationInit = 2.15f;
    float vorticityStrength = 500.f;*/

    // vorticity testing
    /*int numParticles = 30000; 
    float gravity = 0.f; // 4500
    float divergenceModifier = 21.f;
    float gridSize = 74.f; 
    int numPressureIters = 30;
    float diffusionRatio = 0.9f; 
    float flipRatio = 0.80f;
    float seperationInit = 2.75f;
    float vorticityStrength = 2.f;*/

    /*int numParticles = 20000; 
    float gravity = 0.f; // 4500
    float divergenceModifier = 21.f;
    float gridSize = 56.f; 
    int numPressureIters = 30;
    float diffusionRatio = 0.9f; 
    float flipRatio = 0.80f;
    float seperationInit = 3.9f;
    float vorticityStrength = 2.f;*/

    // fire hazard 
    /*int numParticles = 50000; 
    float gravity = 4700.f; 
    float divergenceModifier = 8.5f; 
    float gridSize = 100.f; 
    int numPressureIters = 25; 
    float diffusionRatio = 1.05f; 
    float flipRatio = 0.9f;
    float seperationInit = 1.95f;
    float vorticityStrength = 2.f;*/

    // casual
    /*int32_t numParticles = 15000; 
    float gravity = 5500.f; // 5500
    float divergenceModifier = 8.5f; 
    float gridSize = 70.f; 
    int32_t numPressureIters = 30; // 30 
    float diffusionRatio = 0.85f; 
    float flipRatio = 0.9f;
    float seperationInit = 2.7f; // 2.7
    float vorticityStrength = 400.f;// 350*/

    // grid testing
    /*int numParticles = 1000; 
    float gravity = 0.f; // 5500
    float divergenceModifier = 1.f; 
    float gridSize = 20.f; 
    int numPressureIters = 25; 
    float diffusionRatio = 0.95f; 
    float flipRatio = 0.9f;
    float seperationInit = 4.5f; // 3.5 for 50, 4.5 for 40
    float vorticityStrength = 19.5f;*/

    // low
    int numParticles = 5000; 
    float gravity = 5500.f; // 5500
    float divergenceModifier = 6.5f; 
    float gridSize = 50.f; 
    int numPressureIters = 30; 
    float diffusionRatio = 0.75f; 
    float flipRatio = 0.8f;
    float seperationInit = 3.5f; // 3.5 for 50, 4.5 for 40
    float vorticityStrength = 0.f;

    // low
    /*int numParticles = 7000; 
    float gravity = 5500.f; // 5500
    float divergenceModifier = 6.5f; 
    float gridSize = 50.f; 
    int numPressureIters = 20; 
    float diffusionRatio = 0.95f; 
    float flipRatio = 0.9f;
    float seperationInit = 3.5f; // 3.5 for 50, 4.5 for 40
    float vorticityStrength = 19.5f; */

    // -------------------------------------------

    // for fullscreen
    /*sf::VideoMode desktopMode = sf::VideoMode::getDesktopMode();
    int WIDTH = desktopMode.width;
    int HEIGHT = desktopMode.height;*/

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "FLIP Simulation");

    sf::Font font;
    font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");

    sf::Text text;
    text.setFont(font);
    text.setPosition(10, 10);
    text.setFillColor(sf::Color::White);

    sf::Clock deltaClock;

    window.setFramerateLimit(120); // 120
    //window.setMouseCursorVisible(false);

    int frame = 0;
    int fps = 0;

    float interactionRadius = 100.f;

    float interactionStrength;

    bool leftMouseDown = false;
    bool rightMouseDown = false;

    const uint32_t maxThreads = std::thread::hardware_concurrency();

    const uint32_t numThreads = std::min(static_cast<uint32_t>(1), maxThreads); // 11

    tp::ThreadPool thread_pool(numThreads);

    const float overRelaxation = 1.9f; // 1.9

    Fluid fluid = Fluid(WIDTH, HEIGHT, 1.f * HEIGHT / gridSize, numParticles, gravity, divergenceModifier, diffusionRatio, seperationInit, vorticityStrength, flipRatio, overRelaxation, numPressureIters, thread_pool);  //50

    bool justPressed = false;


    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << flipRatio; 

    std::ostringstream oss2;
    oss2 << std::fixed << std::setprecision(0) << gravity; 

    std::ostringstream oss3;
    oss3 << std::fixed << std::setprecision(1) << divergenceModifier; 

    std::ostringstream oss4;
    oss4 << std::fixed << std::setprecision(1) << vorticityStrength; 

    std::ostringstream oss5;
    oss5 << std::fixed << std::setprecision(0) << numPressureIters; 

    bool forceObjectActive = true; 

    float totalDT = 0;
    float numDT = 0;

    while (window.isOpen())
    {
        sf::Time deltaTime = deltaClock.restart();
        float setDT = 1.f / 240.f;
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
                    if (fluid.getFlipRatio() < 0.99f) {
                        fluid.addToFlipRatio(0.01);
                        oss.str("");  
                        oss.clear();
                        oss << std::fixed << std::setprecision(2) << fluid.getFlipRatio(); 
                    }
                }
                else if (event.key.code == sf::Keyboard::S) {
                    if (fluid.getFlipRatio() > 0.01) {
                        fluid.addToFlipRatio(-0.01);
                        oss.str("");  
                        oss.clear();
                        oss << std::fixed << std::setprecision(2) << fluid.getFlipRatio();
                    }
                }
                if (event.key.code == sf::Keyboard::E) {
                    fluid.addToVorticityStrength(10);
                    oss4.str("");  
                    oss4.clear();
                    oss4 << std::fixed << std::setprecision(1) << fluid.getVorticityStrength(); 
                }
                else if (event.key.code == sf::Keyboard::W) {
                    if (fluid.getVorticityStrength() - 10 >= 0.f) {
                        fluid.addToVorticityStrength(-10);
                        oss4.str("");  
                        oss4.clear();
                        oss4 << std::fixed << std::setprecision(1) << fluid.getVorticityStrength();
                    }
                }
                if (event.key.code == sf::Keyboard::P) {
                    fluid.addToNumPressureIters(1);
                    oss5.str("");  
                    oss5.clear();
                    oss5 << std::fixed << std::setprecision(0) << fluid.getNumPressureIters(); 
                }
                else if (event.key.code == sf::Keyboard::O) {
                    if (fluid.getNumPressureIters() > 0) {
                        fluid.addToNumPressureIters(-1);
                        oss5.str("");  
                        oss5.clear();
                        oss5 << std::fixed << std::setprecision(0) << fluid.getNumPressureIters();
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
                /*else if (event.key.code == sf::Keyboard::Num5) {
                    fluid.setRigidObjectActive(false);
                    fluid.setForceObjectActive(false);
                    fluid.setGeneratorActive(false);
                    fluid.setSolidDrawer(false);
                }*/
                else if (event.key.code == sf::Keyboard::G) {
                    fluid.addToGravity(100);
                    oss2.str("");  
                    oss2.clear();
                    oss2 << std::fixed << std::setprecision(0) << fluid.getGravity();
                }
                else if (event.key.code == sf::Keyboard::N) {
                    fluid.addToGravity(-100);
                    oss2.str("");  
                    oss2.clear();
                    oss2 << std::fixed << std::setprecision(0) << fluid.getGravity();
                }
                else if (event.key.code == sf::Keyboard::C) {
                    fluid.addToDivergenceModifier(0.1);
                    oss3.str("");  
                    oss3.clear();
                    oss3 << std::fixed << std::setprecision(1) << fluid.getDivergenceModifier();
                }
                else if (event.key.code == sf::Keyboard::D) {
                    fluid.addToDivergenceModifier(-0.1);
                    oss3.str("");  
                    oss3.clear();
                    oss3 << std::fixed << std::setprecision(1) << fluid.getDivergenceModifier();
                }
                else if (event.key.code == sf::Keyboard::A) {
                    fluid.setNextRenderPattern();
                }
                else if (event.key.code == sf::Keyboard::F) {
                    fluid.setFireActive(!fluid.getFireActive());
                }
                else if (event.key.code == sf::Keyboard::Q) {
                    //std::cout << "incompressibility time: " << fluid.getTimeForInc() / numDT << ", seperation time: " << fluid.getTimeForSeperation() / numDT << ", grid transfer time: " << fluid.getTimeForTrans() / numDT;
                    window.close();
                }
                else if (event.key.code == sf::Keyboard::Y) {
                    fluid.setStop(!fluid.getStop());
                }
                else if (event.key.code == sf::Keyboard::U) {
                    fluid.setStep(true);
                }
                /*else if (event.key.code == sf::Keyboard::Z) {
                    fluid.setRandomVelField();
                }*/
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

                if (fluid.getGeneratorActive()) {
                    fluid.addToGeneratorRadius(20 * mouseWheelDelta);
                }

                if (fluid.getPencilActive()) {
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

        fluid.update(trueDT, window, leftMouseDown, rightMouseDown, justPressed);

        if (!fluid.getStop()) {
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


        if (fluid.getStop()) {
            text.setFillColor(sf::Color::White);
        }

        text.setPosition(WIDTH - 175, 10);
        text.setString(oss.str());
        window.draw(text);

        text.setPosition(WIDTH - 275, 10);
        text.setString(oss2.str());
        window.draw(text);

        text.setPosition(WIDTH - 375, 10);
        text.setString(oss3.str());
        window.draw(text);

        text.setPosition(WIDTH - 475, 10);
        text.setString(oss4.str());
        window.draw(text);

        text.setPosition(WIDTH - 575, 10);
        text.setString(oss5.str());
        window.draw(text);

        window.display();
 
        justPressed = false;
       
    }

    return 0;
}
