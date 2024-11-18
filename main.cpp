#include <iostream>
#include "fluidSolver.h"
#include "button.h"

using namespace sf;
using namespace std;

FluidSolver fluidSolver(600);

// index sort takes particles as input, sorts them and returns an integer array as output // maybe put into fluid solver?
int* indexSort(Particle* particles, int numFluidParticles)
{
    int* cellIndices = new int[numFluidParticles];
    for (size_t i = 0; i < numFluidParticles; i++)
    {
        // compute cell index with (k, l, m) and the bounding box (find x_min and y_min)
        // increment cellIndices
        // accumulate cellInices
    }
    // sort particles with respect to their index with the help of cellIndices
    return cellIndices;
}

int main()
{
    const int WINDOW_WIDTH = 900;
    const int WINDOW_HEIGHT = 900;

    // states
    bool stopSimulation = false;
    bool showNeighbors = false;

    /* setup window */
    RenderWindow window(VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "FluidSolver", Style::Default);
    window.setFramerateLimit(60);   // prevent too much work for GPU
    window.setPosition(Vector2i(10, 10));

    View view = window.getDefaultView();

    /* load font and prepare text */
    Font font;
    if (!font.loadFromFile("arial.ttf"))
    {
        cout << "font not loaded";
    }
    Text instructions("SHORTCUTS   >>   stop: X | restart: left mouse | zoom: mouse wheel | neighbors : N | (graph : D)", font, 15);
    instructions.setFillColor(Color::Green);
    Text* particleLables = new Text[fluidSolver.numFluidParticles];
    for (size_t i = 0; i < fluidSolver.numFluidParticles; i++)
    {
        particleLables[i] = Text("", font, 9);
        particleLables[i].setFillColor(Color::Green);
    }
    Text cflNumber("", font, 15);
    cflNumber.setFillColor(Color::Green);
    cflNumber.setPosition(Vector2f(0, 15));
    float maxVelocity = 0;

    /* Buttons */
    Button b1(Vector2i(50, 20), Vector2i(10, 50), Color::Green, Text("push", font, 15));

    /* allocate memory for the particle shapes */
    CircleShape* drawingCircles = new CircleShape[fluidSolver.numParticles];
    if (!fluidSolver.particles || !drawingCircles)
    {
        cout << "Memory allocation failed.\n";
    }

    /* initialize all particles */
    fluidSolver.initializeFluidParticles(Vector2f(4, 5));
    fluidSolver.initializeBoundaryParticles();

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
                if (event.mouseButton.button == sf::Mouse::Right)
                {
                    Vector2i mousePos = Mouse::getPosition(window);
                    fluidSolver.initializeFluidParticles(Vector2f((float)mousePos.x / WINDOW_WIDTH * 2 / fluidSolver.H, ((float)WINDOW_HEIGHT - (float)mousePos.y) / WINDOW_WIDTH * 2 / fluidSolver.H));
                    maxVelocity = 0;
                }
                if (event.mouseButton.button == sf::Mouse::Left)
                {
                    if (b1.border.contains(Vector2i(Mouse::getPosition(window))))
                    {
                        cout << "pushed\n";
                    }
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
                break;

            default:
                break;
            }
        }


        if (!stopSimulation)
        {
            /* Update (SPH Fluid Solver) */
            fluidSolver.neighborSearchNN(2);
            fluidSolver.computeDensityAndPressure();
            fluidSolver.computeAccelerations();
            fluidSolver.updatePositions();
        }

        /* Draw */

        window.setView(view);
        window.clear(); // don't draw on top of the previous frame

        for (size_t i = 0; i < fluidSolver.numParticles; i++)
        {
            drawingCircles[i].setRadius(fluidSolver.H / 2.f * WINDOW_WIDTH / 2.f);    // h is defined as the "diameter"
            drawingCircles[i].setPosition(Vector2f((fluidSolver.particles[i].position.x + 1.f) * WINDOW_WIDTH / 2.f, WINDOW_HEIGHT - (fluidSolver.particles[i].position.y + 1.f) * WINDOW_WIDTH / 2.f));   // the shapes to be drawn have to be updated independently, scale
            if (i < fluidSolver.numFluidParticles)
            {
                drawingCircles[i].setFillColor(Color::Blue);
            }
            window.draw(drawingCircles[i]);
        }

        /* text */
        if (showNeighbors)
        {
            for (size_t i = 0; i < fluidSolver.numFluidParticles; i++)
            {
                particleLables[i].setString(to_string(fluidSolver.particles[i].neighbors.size()));
                Vector2f pixelCoord = Vector2f((fluidSolver.particles[i].position.x + 1.f) * WINDOW_WIDTH / 2.f, WINDOW_HEIGHT - (fluidSolver.particles[i].position.y + 1.f) * WINDOW_WIDTH / 2.f);
                particleLables[i].setPosition(pixelCoord);
                window.draw(particleLables[i]);
            }
        }

        for (size_t i = 0; i < fluidSolver.numFluidParticles; i++)
        {
            maxVelocity = max(maxVelocity, sqrt(fluidSolver.particles[i].velocity.x * fluidSolver.particles[i].velocity.x + fluidSolver.particles[i].velocity.y * fluidSolver.particles[i].velocity.y));
        }

        cflNumber.setString("CFL: lambda >= " + to_string((fluidSolver.TIME_STEP * maxVelocity) / fluidSolver.H) + ", maxTimeStep: " + to_string(fluidSolver.H / maxVelocity));
        window.draw(cflNumber);
        window.draw(instructions);

        /* buttons */
        window.draw(b1.shape);
        window.draw(b1.name);

        /* Display */
        window.display();
    }

    /* deallocate memory */
    delete[] drawingCircles;
    delete[] particleLables;

    return EXIT_SUCCESS;
}