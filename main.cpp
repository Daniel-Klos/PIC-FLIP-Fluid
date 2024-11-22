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
    // Adjust the sim to a good size for you! Note that if you change these numbers, you will have to change the settings!
    int WIDTH = 2500; 
    int HEIGHT = 1300; 

    // SETTINGS
    // ---------------------------------------

    // 1) space
    /*int numParticles = 1000; 
    float gravity = 0.f; 
    float divergenceModifier = 15.f; 
    float gridSize = 90.f; 
    int numPressureIters = 30; 
    float diffusionRatio = 0.8f;
    float flipRatio = 0.8f;*/

    // 2) casual
    int numParticles = 15000; 
    float gravity = 5500.f; 
    float divergenceModifier = 8.5f; 
    float gridSize = 70.f; 
    int numPressureIters = 25; 
    float diffusionRatio = 0.95f; 
    float flipRatio = 0.9f;
    float seperationInit = 2.7f;
    float vorticityStrength = 2.f; // messes with the grid in some irreversible way at the start of the sim if set to anything > 0

    // 3) lots
    /*int numParticles = 25000; 
    float gravity = 4700.f; 
    float divergenceModifier = 15.5f; 
    float gridSize = 75.f; 
    int numPressureIters = 25; 
    float diffusionRatio = 0.95f; 
    float flipRatio = 0.9f;*/

    // 4) a lot
    /*int numParticles = 30000; 
    float gravity = 4500.f; 
    float divergenceModifier = 
    float gridSize = 90.f; 
    int numPressureIters = 30;
    float diffusionRatio = 0.9f; 
    float flipRatio = 0.80f;*/

    // fire hazard (set dt to a fixed amount for this aint no way this is running in real time)
    /*int numParticles = 50000; 
    float gravity = 4700.f; 
    float divergenceModifier = 8.5f; 
    float gridSize = 120.f; 
    int numPressureIters = 25; 
    float diffusionRatio = 1.05f; 
    float flipRatio = 0.9f;
    float seperationInit = 1.55f;*/

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

    //window.setFramerateLimit(120);
    //window.setMouseCursorVisible(false);

    int frame = 0;
    int fps = 0;

    float interactionRadius = 100.f;

    float interactionStrength;

    bool leftMouseDown = false;
    bool rightMouseDown = false;

    const uint32_t maxThreads = std::thread::hardware_concurrency();

    const uint32_t numThreads = maxThreads < 11 ? maxThreads : 11;

    tp::ThreadPool thread_pool(numThreads);

    Fluid fluid = Fluid(WIDTH, HEIGHT, 1.f * HEIGHT / gridSize, numParticles, gravity, divergenceModifier, diffusionRatio, seperationInit, vorticityStrength, thread_pool);  //50

    const float overRelaxation = 1.91f;

    bool justPressed = false;


    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << flipRatio; 

    std::ostringstream oss2;
    oss2 << std::fixed << std::setprecision(0) << gravity; 

    std::ostringstream oss3;
    oss3 << std::fixed << std::setprecision(1) << divergenceModifier; 

    std::ostringstream oss4;
    oss4 << std::fixed << std::setprecision(1) << vorticityStrength; 

    bool forceObjectActive = true; 

    float totalDT = 0;
    float numDT = 0;

    while (window.isOpen())
    {
        sf::Time deltaTime = deltaClock.restart();
        float dt = deltaTime.asSeconds();

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
                    if (flipRatio < 0.99) {
                        flipRatio += 0.01;
                        oss.str("");  
                        oss.clear();
                        oss << std::fixed << std::setprecision(2) << flipRatio; 
                    }
                }
                else if (event.key.code == sf::Keyboard::S) {
                    if (flipRatio > 0.01) {
                        flipRatio -= 0.01;
                        oss.str("");  
                        oss.clear();
                        oss << std::fixed << std::setprecision(2) << flipRatio;
                    }
                }
                if (event.key.code == sf::Keyboard::E) {
                    vorticityStrength += 0.5;
                    oss4.str("");  
                    oss4.clear();
                    oss4 << std::fixed << std::setprecision(1) << vorticityStrength; 
                }
                else if (event.key.code == sf::Keyboard::W) {
                    if (vorticityStrength - 0.5 > 0.f) {
                        vorticityStrength -= 0.5;
                        oss4.str("");  
                        oss4.clear();
                        oss4 << std::fixed << std::setprecision(1) << vorticityStrength;
                    }
                }
                else if (event.key.code == sf::Keyboard::Num1) {
                    fluid.setRigidObjectActive(true);
                    fluid.setForceObjectActive(false);
                    fluid.setGeneratorActive(false);
                }
                else if (event.key.code == sf::Keyboard::Num2) {
                    fluid.setRigidObjectActive(false);
                    fluid.setForceObjectActive(true);
                    fluid.setGeneratorActive(false);
                }
                else if (event.key.code == sf::Keyboard::Num3) {
                    fluid.setRigidObjectActive(false);
                    fluid.setForceObjectActive(false);
                    fluid.setGeneratorActive(true);
                }
                else if (event.key.code == sf::Keyboard::T) {
                    fluid.addToForceObjectRadius(10);
                    fluid.addToGeneratorRadius(10);
                }
                else if (event.key.code == sf::Keyboard::R) {
                    fluid.addToForceObjectRadius(-10);
                    fluid.addToGeneratorRadius(-10);
                }
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
                else if (event.key.code == sf::Keyboard::F) {
                    fluid.setDiffuseColors(!fluid.getDiffuseColors());
                }
                else if (event.key.code == sf::Keyboard::Q) {
                    //std::cout << "incompressibility time: " << fluid.getTimeForInc() / numDT << ", seperation time: " << fluid.getTimeForSeperation() / numDT << ", grid transfer time: " << fluid.getTimeForTrans() / numDT;
                    window.close();
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
        }

        window.clear();

        fluid.simulate(dt, window, leftMouseDown, justPressed, numPressureIters, overRelaxation, flipRatio, rightMouseDown);

        frame++;
        if (frame == 30) {
            fps = (int)(1.f / dt);
            frame = 0;
        }

        text.setPosition(WIDTH - 70, 10);
        text.setString(std::to_string(fps)); 
        window.draw(text);


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

        window.display();
 
        justPressed = false;
       
    }

    return 0;
}
