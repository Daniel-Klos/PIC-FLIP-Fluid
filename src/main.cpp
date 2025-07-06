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
    int numParticles = 30000; // 20000 --- start with a kinda large number so that good density sample is taken at the start of the simulation
    float gravityY = 5500.f; // 5500
    float gravityX = 0.f;
    int gridNumX = 348; // 350
    float diffusionRatio = 0.75f; 
    float flipRatio = 0.9f; // 0.9f
    float vorticityStrength = 0.f;
    int maxFps = 120;
    uint32_t numThreads = 6; // 6

    const uint32_t maxThreads = std::thread::hardware_concurrency();

    numThreads = std::min(numThreads, maxThreads); // 10, 11, 16

    tp::ThreadPool thread_pool(numThreads);

    FrameContext frame_context = FrameContext(WIDTH, HEIGHT);

    FluidState fluid_attributes = FluidState(numParticles, gridNumX, vorticityStrength, flipRatio, gravityX, gravityY, frame_context, thread_pool);

    sf::RenderWindow window(sf::VideoMode(fluid_attributes.frame_context.WIDTH, fluid_attributes.frame_context.HEIGHT), "FLIP Simulation");

    SceneHandler scene_handler = SceneHandler(fluid_attributes, window, maxFps);

    while(window.isOpen()) {
        scene_handler.simulate();
    }

    return 0;
}