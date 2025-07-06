#pragma once
#include <vector>
#include <numeric>
#include <iostream>
#include <SFML/Graphics.hpp>
#include "simulation_state.hpp"
#include "utils.hpp"

class SceneRenderer;

struct ObstacleRenderer {
    FluidState &fluid_attributes;

    sf::VertexArray obstacleVa{sf::PrimitiveType::Quads};

    sf::Texture obstacleTexture;
    sf::Vector2f obstacle_texture_size;

    sf::RenderStates obstacleStates;

    sf::RenderWindow &window;

    sf::RectangleShape pencil;

    float pencilSeparationX;
    float pencilSeparationY;
    int n;

    ObstacleRenderer(FluidState &fas, sf::RenderWindow &w): fluid_attributes(fas), window(w) {
        pencil.setSize(sf::Vector2f{fluid_attributes.cellSpacing, fluid_attributes.cellSpacing});
        pencil.setOrigin(fluid_attributes.halfSpacing, fluid_attributes.halfSpacing);
        pencil.setOutlineThickness(1);
        pencil.setOutlineColor(sf::Color::Black);
        pencilSeparationX = fluid_attributes.cellSpacing;
        pencilSeparationY = fluid_attributes.cellSpacing;
    }

    void render_obstacles() {
        window.draw(obstacleVa, obstacleStates);
    }

    void drawPencil(int pencilRadius) {
        sf::Vector2i localPos = sf::Vector2i{fluid_attributes.frame_context.world_mouse_pos / fluid_attributes.cellSpacing};

        sf::Vector2f simPos = sf::Vector2f(localPos) * fluid_attributes.cellSpacing + sf::Vector2f{fluid_attributes.halfSpacing, fluid_attributes.halfSpacing};

        sf::Vector2f screenPos = (simPos - fluid_attributes.frame_context.offset) * fluid_attributes.frame_context.zoom_amount + fluid_attributes.frame_context.center;

        if (fluid_attributes.frame_context.leftMouseDown || !fluid_attributes.frame_context.rightMouseDown) {
            pencil.setFillColor(sf::Color(0, 150, 0));
        }
        else {
            pencil.setFillColor(sf::Color(150, 0, 0));
        }

        for (int i = -pencilRadius; i <= pencilRadius; ++i) {
            for (int j = -pencilRadius; j <= pencilRadius; ++j) {
                pencil.setPosition(screenPos.x + i * pencilSeparationX, screenPos.y + j * pencilSeparationY);
                window.draw(pencil);
            }
        }
    }

};