#pragma once
#include "simulation_state.hpp"
#include "fluid_rendering.hpp"
#include "obstacle_rendering.hpp"

struct SceneRenderer {
    FluidState &fluid_attributes;
    FluidRenderer fluid_renderer;
    ObstacleRenderer obstacle_renderer;
    sf::RenderWindow &window;

    SceneRenderer(FluidState &fas, sf::RenderWindow& w): fluid_attributes(fas), fluid_renderer(fas, w), obstacle_renderer(fas, w), window(w) {}

    void render_fluid() {
        if (fluid_renderer.renderPattern == 0) {
            fluid_renderer.UpdateVaDiffusionMulti();
        }

        else if (fluid_renderer.renderPattern == 1) {
            fluid_renderer.UpdateVaVelocityMulti();
        }

        else if (fluid_renderer.renderPattern == 2) {
            fluid_renderer.UpdateVaVorticityMulti();
        }

        else if (fluid_renderer.renderPattern == 3) {
            fluid_renderer.UpdateVaTemperatureMulti();
        }
        else if (fluid_renderer.renderPattern == 4) {
            fluid_renderer.UpdateDivergenceVaMulti();
            fluid_renderer.DrawDivergences();
            //fluid_renderer.drawActiveUVNodes(window);
        }
        else if (fluid_renderer.renderPattern == 5) {
            fluid_renderer.UpdateVaCustomMulti();
        }

        if (fluid_renderer.renderPattern != 4) {
            fluid_renderer.DrawParticles();
        }
    }

};