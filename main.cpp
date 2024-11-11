#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <array>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <thread>

#include "fluid.hpp"

int main()
{
    // make gravity larger for faster seeming simulations, for 12000 particles like 6500 gravity and 9.3 divergence modifier
    int numParticles = 15000; 
    float gravity = 4500.f; // 4500
    // gs 75, pc 25k : 11.3, gs 75, pc 30k : 20, gs 75, pc 20k : 5
    float divergenceModifier = 8.f; // 13 for a wide simulation: 2.7 for 8k, ~3.5 for 10k, 3.9 for 12k, 4.7-5-7.6 for 15k  |  for tall: 15 for 10k   (also depends on grid size)
    float gridSize = 55; // for wide: 50-70 (65) (75)  |  for tall: 80-90
    int numPressureIters = 25; // for a wide simulation : ~20  |  for a tall simulation: ~50
    float diffusionRatio = 0.85f; // 1.05

    // for fullscreen
    /*sf::VideoMode desktopMode = sf::VideoMode::getDesktopMode();
    int WIDTH = desktopMode.width;
    int HEIGHT = desktopMode.height;*/

    int WIDTH = 2500; //wide simulation: 2000  |  tall simulation: 800
    int HEIGHT = 1000; // wide simulatoin: 1000  |  tall simulation: 1300 

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

    const uint32_t numThreads = maxThreads < 10 ? maxThreads : 10;

    tp::ThreadPool thread_pool(numThreads);

    Fluid fluid = Fluid(1000, WIDTH, HEIGHT, 1.f * HEIGHT / gridSize, numParticles, gravity, divergenceModifier, diffusionRatio, thread_pool);  //50

    float overRelaxation = 1.91f;
    float flipRatio = 0.90f;

    bool justPressed = false;


    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << flipRatio; 

    std::ostringstream oss2;
    oss2 << std::fixed << std::setprecision(0) << gravity; 

    std::ostringstream oss3;
    oss3 << std::fixed << std::setprecision(1) << divergenceModifier; 

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
                else if (event.key.code == sf::Keyboard::Num1) {
                    forceObjectActive = false;
                }
                else if (event.key.code == sf::Keyboard::Num2) {
                    forceObjectActive = true;
                }
                else if (event.key.code == sf::Keyboard::T) {
                    fluid.addToForceObjectRadius(3);
                }
                else if (event.key.code == sf::Keyboard::R) {
                    fluid.addToForceObjectRadius(-3);
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

        fluid.simulate(dt, window, forceObjectActive, leftMouseDown, justPressed, numPressureIters, overRelaxation, flipRatio, rightMouseDown);

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

        window.display();
 
        justPressed = false;
       
    }

    return 0;
}
