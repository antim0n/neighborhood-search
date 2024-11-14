#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>

using namespace sf;

class Button
{
public:
	RectangleShape shape;
	IntRect border;
	Text name;
	Button(Vector2i size, Vector2i position, Color color, Text text);
};