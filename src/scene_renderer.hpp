#pragma once
#include "simulation_state.hpp"
#include "fluid_rendering.hpp"
#include "obstacle_rendering.hpp"

struct SceneRenderer {
    FluidState &fluid_attributes;
    FluidRenderer fluid_renderer;
    ObstacleRenderer obstacle_renderer;
    sf::RenderWindow &window;

    bool zoomObjectActive = false;
    float zoom_amount = 1.f;
    sf::Vector2f prev_world_mouse_pos;
    sf::Vector2f offset;
    sf::Vector2f center;
    sf::Transform transform;

    SceneRenderer(FluidState &fas, sf::RenderWindow& w): fluid_attributes(fas), fluid_renderer(fas, w), obstacle_renderer(fas, w), window(w) {
        float centerX = fluid_attributes.frame_context.WIDTH / 2.f;
        float centerY = fluid_attributes.frame_context.HEIGHT / 2.f;
        offset = sf::Vector2f{centerX, centerY};
        center = sf::Vector2f{centerX, centerY};
    }

    // pass in sf::Vector
    template <typename T>
    T screenToWorld(T screen_pos) const {
        return offset + (screen_pos - center) / zoom_amount;
    }

    void render_scene() {
        fluid_renderer.render_fluid();
        obstacle_renderer.render_obstacles();
    }

    void zoomInto(sf::Vector2f zoom_position, float f) {
        fluid_attributes.frame_context.world_mouse_pos = screenToWorld(zoom_position);

        zoom_amount *= f;

        offset = fluid_attributes.frame_context.world_mouse_pos - (zoom_position - center) / zoom_amount;

        zoom_scene();
    }

     void zoom_scene() {
        const float z = zoom_amount;
        transform = sf::Transform::Identity;
        transform.translate(center);
        transform.scale(z, z);
        transform.translate(-offset);

        fluid_renderer.fluidStates.transform = transform;
        fluid_renderer.cellStates.transform = transform;
        
        obstacle_renderer.obstacleStates.transform = transform;

        normalize_objects();
    }

    void reset_zoom() {
        zoom_amount = 1.f;
        offset = center;
        fluid_attributes.frame_context.world_mouse_pos = center;
        zoom_scene();
    }
    
    void wheelZoom(int delta) {
        const float zoom_mag = 1.09f;
        const float delta_new = delta > 0 ? zoom_mag : 1.0f / zoom_mag;
        zoomInto(fluid_attributes.frame_context.screen_mouse_pos, delta_new);
    }

    void dragCamera() {

    }

    void normalize_objects() {
        obstacle_renderer.pencil.setSize(sf::Vector2f{fluid_attributes.cellSpacing * zoom_amount, fluid_attributes.cellSpacing * zoom_amount});
        obstacle_renderer.pencil.setOrigin(fluid_attributes.halfSpacing * zoom_amount, fluid_attributes.halfSpacing * zoom_amount);
        obstacle_renderer.pencilSeparationX = std::max(fluid_attributes.cellSpacing * zoom_amount, 0.1f);
        obstacle_renderer.pencilSeparationY = std::max(fluid_attributes.cellSpacing * zoom_amount, 0.1f);
    }

    bool getZoomObjectActive() {
        return zoomObjectActive;
    }

    void setZoomObjectActive(bool set) {
        zoomObjectActive = set;
    }

};