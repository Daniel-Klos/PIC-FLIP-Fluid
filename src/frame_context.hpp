#pragma once

struct FrameContext {
    float dt;
    bool leftMouseDown;
    bool rightMouseDown;
    bool justPressed;
    float WIDTH;
    float HEIGHT;

    sf::Vector2f screen_mouse_pos;
    sf::Vector2f world_mouse_pos;
    sf::Vector2f simulation_mouse_pos;

    FrameContext(float WIDTH_, float HEIGHT_): WIDTH(WIDTH_), HEIGHT(HEIGHT_) {}
};