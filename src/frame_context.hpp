#pragma once

struct FrameContext {
    float dt;
    int maxFps;
    
    bool leftMouseDown = false;
    bool rightMouseDown = false;
    bool justPressed = false;
    float WIDTH;
    float HEIGHT;

    sf::Vector2f screen_mouse_pos;
    sf::Vector2f world_mouse_pos;

    sf::Vector2f prev_screen_mouse_pos;
    sf::Vector2f prev_world_mouse_pos;

    sf::Vector2f offset;
    sf::Vector2f center;
    sf::Transform transform;
    float zoom_amount = 1.f;
    bool zooming = false;

    FrameContext(float WIDTH_, float HEIGHT_): WIDTH(WIDTH_), HEIGHT(HEIGHT_) {}
};