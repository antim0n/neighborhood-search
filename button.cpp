#include "button.h"

Button::Button(Vector2i size, Vector2i position, Color color, Text text)
{
	shape = RectangleShape(Vector2f(size.x, size.y));
	shape.setPosition(Vector2f(position.x, position.y));
	shape.setFillColor(color);
	name = text;
	name.setPosition(Vector2f(position.x + size.x / 10, position.y));
	border = IntRect(position, size);
}