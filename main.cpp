#include <iostream>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>

#include "fluidSolver.h"

using namespace sf;
using namespace std;

const int NUM_FLUID_PARTICLES = 600;
const int NUM_BOUNDARY_PARTICLES = 460;
const int NUMBER_OF_PARTICLES = NUM_FLUID_PARTICLES + NUM_BOUNDARY_PARTICLES;

const int WINDOW_WIDTH = 900;
const int WINDOW_HEIGHT = 900;

// states
bool stopSimulation = false;
bool showNeighbors = false;

static Vector2f particleCoordsToPixel(Vector2f position)
{
    return Vector2f((position.x + 1.f) * WINDOW_WIDTH / 2.f, WINDOW_HEIGHT - (position.y + 1.f) * WINDOW_WIDTH / 2.f);
}


int main()
{
    /* setup window */
    RenderWindow window(VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "FluidSolver", Style::Default);
    window.setFramerateLimit(60);   // prevent too much work for GPU
    window.setPosition(Vector2i(10, 10));

    View view = window.getDefaultView();

    /* load font */
    Font font;
    if (!font.loadFromFile("arial.ttf"))
    {
        cout << "font not loaded";
    }
    Text text;
    text.setFont(font);
    text.setString("SHORTCUTS   >>   stop: X | restart: left mouse | zoom: mouse wheel | neighbors : N | (graph : D)");
    text.setCharacterSize(15);
    text.setFillColor(Color::Green);

    Text* particleLables = new Text[NUM_FLUID_PARTICLES]; // TODO: create extra text setup function?
    for (size_t i = 0; i < NUM_FLUID_PARTICLES; i++)
    {
        particleLables[i].setFont(font);
        particleLables[i].setCharacterSize(9);
        particleLables[i].setFillColor(Color::Green);
    }
    Text cflNumber;
    float maxVelocity = 0;
    cflNumber.setFont(font);
    cflNumber.setCharacterSize(15);
    cflNumber.setFillColor(Color::Green);
    cflNumber.setPosition(Vector2f(0, 15));

    /* allocate memory for the particles and their shapes */
    Particle* particles = new Particle[NUMBER_OF_PARTICLES];
    CircleShape* drawingCircles = new CircleShape[NUMBER_OF_PARTICLES];

    if (!particles || !drawingCircles)
    {
        cout << "Memory allocation failed.\n";
    }

    /* initialize all fluid particles */
    initializeFluidParticles(particles, NUM_FLUID_PARTICLES, Vector2f(4, 5));

    /* initialize all boundary particles */
    initializeBoundaryParticles(particles, NUM_FLUID_PARTICLES, NUMBER_OF_PARTICLES);

    /* simulation and rendering loop */
    while (window.isOpen())
    {
        Event event;
        while (window.pollEvent(event)) // check for user inputs (only for closing the window)
        {
            switch (event.type)
            {
            case Event::Closed:
                window.close();
                break;

            case Event::KeyPressed:
                if (event.key.scancode == sf::Keyboard::Scan::X)
                {
                    stopSimulation = !stopSimulation;
                }
                else if (event.key.scancode == sf::Keyboard::Scan::N)
                {
                    showNeighbors = !showNeighbors;
                }
                else if (event.key.scancode == sf::Keyboard::Scan::D)
                {
                    /*drawGraphs = !drawGraphs;
                    myTime = 0;*/
                }
                break;

            case Event::MouseButtonPressed:
                if (event.mouseButton.button == sf::Mouse::Left)
                {
                    Vector2i mousePos = Mouse::getPosition(window);
                    initializeFluidParticles(particles, NUM_FLUID_PARTICLES, Vector2f((float)mousePos.x / WINDOW_WIDTH * 2 / H, ((float)WINDOW_HEIGHT - (float)mousePos.y) / WINDOW_WIDTH * 2 / H));
                    maxVelocity = 0;
                }
                break;

            case Event::MouseWheelScrolled:
                if (event.mouseWheelScroll.delta >= 1)
                {
                    view.zoom(0.95f); // TODO: fitted scaling
                }
                if (event.mouseWheelScroll.delta <= -1)
                {
                    view.zoom(1.05f);
                }

            default:
                break;
            }
        }


        if (!stopSimulation)
        {
            /* Update (SPH Fluid Solver) */
            neighborSearchNN(particles, NUM_FLUID_PARTICLES, NUMBER_OF_PARTICLES, 2);
            computeDensityAndPressure(particles, NUM_FLUID_PARTICLES);
            computeAccelerations(particles, NUM_FLUID_PARTICLES);
            updatePositions(particles, NUM_FLUID_PARTICLES);
        }

        /* Draw */

        window.setView(view);
        window.clear(); // don't draw on top of the previous frame

        for (size_t i = 0; i < NUMBER_OF_PARTICLES; i++)
        {
            drawingCircles[i].setRadius(H / 2.f * WINDOW_WIDTH / 2.f);    // h is defined as the "diameter"
            drawingCircles[i].setPosition(Vector2f((particles[i].position.x + 1.f) * WINDOW_WIDTH / 2.f, WINDOW_HEIGHT - (particles[i].position.y + 1.f) * WINDOW_WIDTH / 2.f));   // the shapes to be drawn have to be updated independently, scale
            if (i < NUM_FLUID_PARTICLES)
            {
                drawingCircles[i].setFillColor(Color::Blue);
            }
            window.draw(drawingCircles[i]);
        }

        /* text */
        if (showNeighbors)
        {
            for (size_t i = 0; i < NUM_FLUID_PARTICLES; i++)
            {
                particleLables[i].setString(to_string(particles[i].neighbors.size()));
                particleLables[i].setPosition(particleCoordsToPixel(particles[i].position));
                window.draw(particleLables[i]);
            }
        }

        for (size_t i = 0; i < NUM_FLUID_PARTICLES; i++)
        {
            maxVelocity = max(maxVelocity, sqrt(particles[i].velocity.x * particles[i].velocity.x + particles[i].velocity.y * particles[i].velocity.y));
        }

        cflNumber.setString("CFL: lambda >= " + to_string((TIME_STEP * maxVelocity) / H) + ", maxTimeStep: " + to_string(H / maxVelocity));
        window.draw(cflNumber);
        window.draw(text);

        /* Display */
        window.display();
    }

    /* deallocate memory */
    delete[] particles;
    delete[] drawingCircles;
    delete[] particleLables;

    return EXIT_SUCCESS;
}