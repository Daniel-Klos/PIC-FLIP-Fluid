#pragma once

struct FrameContext {
    float dt;
    float mouseX;
    float mouseY;
    bool leftMouseDown;
    bool rightMouseDown;
    bool justPressed;
    float WIDTH;
    float HEIGHT;

    FrameContext(float WIDTH_, float HEIGHT_): WIDTH(WIDTH_), HEIGHT(HEIGHT_) {}
};