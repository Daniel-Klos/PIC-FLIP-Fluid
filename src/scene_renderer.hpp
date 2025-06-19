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

    SceneRenderer(FluidState &fas, sf::RenderWindow& w): fluid_attributes(fas), fluid_renderer(fas, w), obstacle_renderer(fas, w), window(w) {
        float centerX = fluid_attributes.frame_context.WIDTH / 2.f;
        float centerY = fluid_attributes.frame_context.HEIGHT / 2.f;
        fluid_attributes.frame_context.offset = sf::Vector2f{centerX, centerY};
        fluid_attributes.frame_context.center = sf::Vector2f{centerX, centerY};

        fluid_attributes.frame_context.zoom_amount = 1.f;
    }

    // pass in sf::Vector
    template <typename T>
    T screenToWorld(T screen_pos) const {
        return fluid_attributes.frame_context.offset + (screen_pos - fluid_attributes.frame_context.center) / fluid_attributes.frame_context.zoom_amount;
    }

    void render_scene() {
        fluid_renderer.render_fluid();
        obstacle_renderer.render_obstacles();
    }

    void zoomInto(sf::Vector2f zoom_position, float f) {
        fluid_attributes.frame_context.world_mouse_pos = screenToWorld(zoom_position);

        fluid_attributes.frame_context.zoom_amount *= f;

        fluid_attributes.frame_context.offset = fluid_attributes.frame_context.world_mouse_pos - (zoom_position - fluid_attributes.frame_context.center) / fluid_attributes.frame_context.zoom_amount;

        zoom_scene();
    }

     void zoom_scene() {
        fluid_attributes.frame_context.zooming = true;
        const float z = fluid_attributes.frame_context.zoom_amount;
        fluid_attributes.frame_context.transform = sf::Transform::Identity;
        fluid_attributes.frame_context.transform.translate(fluid_attributes.frame_context.center);
        fluid_attributes.frame_context.transform.scale(z, z);
        fluid_attributes.frame_context.transform.translate(-fluid_attributes.frame_context.offset);

        fluid_renderer.fluidStates.transform = fluid_attributes.frame_context.transform;
        fluid_renderer.cellStates.transform = fluid_attributes.frame_context.transform;
        
        obstacle_renderer.obstacleStates.transform = fluid_attributes.frame_context.transform;

        normalize_objects();
    }

    void reset_zoom() {
        fluid_attributes.frame_context.zoom_amount = 1.f;
        fluid_attributes.frame_context.offset = fluid_attributes.frame_context.center;
        fluid_attributes.frame_context.world_mouse_pos = fluid_attributes.frame_context.center;
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
        obstacle_renderer.pencil.setSize(sf::Vector2f{fluid_attributes.cellSpacing * fluid_attributes.frame_context.zoom_amount, fluid_attributes.cellSpacing * fluid_attributes.frame_context.zoom_amount});
        obstacle_renderer.pencil.setOrigin(fluid_attributes.halfSpacing * fluid_attributes.frame_context.zoom_amount, fluid_attributes.halfSpacing * fluid_attributes.frame_context.zoom_amount);
        obstacle_renderer.pencilSeparationX = std::max(fluid_attributes.cellSpacing * fluid_attributes.frame_context.zoom_amount, 0.1f);
        obstacle_renderer.pencilSeparationY = std::max(fluid_attributes.cellSpacing * fluid_attributes.frame_context.zoom_amount, 0.1f);
    }

    bool getZoomObjectActive() {
        return zoomObjectActive;
    }

    void setZoomObjectActive(bool set) {
        zoomObjectActive = set;
    }

};