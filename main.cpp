#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include <iostream>
#include <chrono>
#include "fluidSolver.h"
#include "neighborSearch.h"
#include "button.h"

using namespace std;
using namespace std::chrono;
using namespace sf;

FluidSolver fluidSolver(100000);

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

    // fps calculation
    Clock clock;
    Time previousTime = clock.getElapsedTime();
    Time currentTime;
    int fps = 60;
    int counter = 0;

    /* load font and prepare text */
    Font font;
    if (!font.loadFromFile("arial.ttf"))
    {
        cout << "font not loaded";
    }
    Text instructions("SHORTCUTS   >>   stop: X | restart: left mouse | zoom: mouse wheel | neighbors : N | (graph : D)", font, 15);
    instructions.setFillColor(Color::Green);
    Text* particleLables = new Text[fluidSolver.numParticles];
    for (size_t i = 0; i < fluidSolver.numParticles; i++)
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
    fluidSolver.initializeFluidParticles(Vector2f(4.f * fluidSolver.H, 5.f * fluidSolver.H));
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
                    Vector2f newOffset = Vector2f(static_cast<float>(mousePos.x) / (static_cast<float>(WINDOW_HEIGHT) / 100.f),
                        static_cast<float>(WINDOW_HEIGHT - mousePos.y) / (static_cast<float>(WINDOW_HEIGHT) / 100.f)); // pixel to particle
                    fluidSolver.initializeFluidParticles(newOffset);
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
            // runtime measurement
            auto start = high_resolution_clock::now();

            indexSortConstruction(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            indexSortQuery(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);

            /*for (size_t i = 0; i < fluidSolver.numFluidParticles; i++)
            {
                for (size_t j = 0; j < fluidSolver.particles[i].neighbors.size(); j++)
                {
                    cout << fluidSolver.particles[i].neighbors[j]->index << " ";
                }
                cout << endl;
            }
            cout << endl;*/

            // fluidSolver.neighborSearchNN(2);

            /*for (size_t i = 0; i < fluidSolver.numFluidParticles; i++)
            {
                for (size_t j = 0; j < fluidSolver.particles[i].neighbors.size(); j++)
                {
                    cout << fluidSolver.particles[i].neighbors[j]->index << " ";
                }
                cout << endl;
            }*/

            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<chrono::milliseconds>(stop - start);
            cout << duration.count() << endl;
            cout << chrono::duration<double>(duration).count() << endl;

            /*fluidSolver.computeDensityAndPressure();
            fluidSolver.computeAccelerations();
            fluidSolver.updatePositions();*/
            // stopSimulation = true;
        }

        /* Draw */

        window.setView(view);
        window.clear(); // don't draw on top of the previous frame

        for (size_t i = 0; i < fluidSolver.numParticles; i++)
        {
            drawingCircles[i].setRadius((fluidSolver.H / 2.f) * (static_cast<float>(WINDOW_HEIGHT) / 100.f));    // h is defined as the "diameter"
            drawingCircles[i].setPosition(Vector2f((fluidSolver.particles[i].position.x) * (static_cast<float>(WINDOW_WIDTH) / 100.f),
                static_cast<float>(WINDOW_HEIGHT) - fluidSolver.particles[i].position.y * (static_cast<float>(WINDOW_WIDTH) / 100.f)));   // the shapes to be drawn have to be updated independently, scale
            if (fluidSolver.particles[i].isFluid)
            {
                drawingCircles[i].setFillColor(Color::Blue);
            }
            else
            {
                drawingCircles[i].setFillColor(Color::White);
            }
            window.draw(drawingCircles[i]);
        }

        // grid
        for (size_t i = 0; i < boundingBox[4] + 1; i++) // vertical
        {
            Vertex line[] = {
                Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[0] + 2.f * fluidSolver.H * i, boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Yellow),
                Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[0] + 2.f * fluidSolver.H * i, boundingBox[2] + 2.f * fluidSolver.H * (boundingBox[5])), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Yellow) };
            window.draw(line, 2, Lines);
        }
        for (size_t i = 0; i < boundingBox[5] + 1; i++) // horizontal
        {
            Vertex line[] = {
                Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[2] + 2.f * fluidSolver.H * i), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Yellow),
                Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[0] + 2.f * fluidSolver.H * boundingBox[4], boundingBox[2] + 2.f * fluidSolver.H * i), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Yellow) };
            window.draw(line, 2, Lines);
        }

        // bounding box
        Vertex line1[] = {
            Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Red),
            Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[1], boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Red) };
        Vertex line2[] = {
            Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Red),
            Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[3]), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Red) };
        Vertex line3[] = {
            Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[3]), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Red),
            Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[1], boundingBox[3]), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Red) };
        Vertex line4[] = {
            Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[1], boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Red),
            Vertex(fluidSolver.particleToPixelCoord(Vector2f(boundingBox[1], boundingBox[3]), WINDOW_WIDTH, WINDOW_HEIGHT), Color::Red) };
        window.draw(line1, 2, Lines);
        window.draw(line2, 2, Lines);
        window.draw(line3, 2, Lines);
        window.draw(line4, 2, Lines);

        /* text */
        if (showNeighbors)
        {
            for (size_t i = 0; i < fluidSolver.numParticles; i++)
            {
                if (fluidSolver.particles[i].isFluid)
                {
                    particleLables[i].setString(to_string(fluidSolver.particles[i].neighbors.size()));
                    Vector2f pixelCoord = Vector2f((fluidSolver.particles[i].position.x) * (static_cast<float>(WINDOW_WIDTH) / 100.f),
                        static_cast<float>(WINDOW_HEIGHT) - fluidSolver.particles[i].position.y * (static_cast<float>(WINDOW_WIDTH) / 100.f)); // particle to pixel
                    particleLables[i].setPosition(pixelCoord);
                    window.draw(particleLables[i]);
                }
            }
        }

        // fps calculation
        currentTime = clock.getElapsedTime();
        counter += 1;
        if (counter > 20)
        {
            fps = static_cast<int>(1.0f / (currentTime.asSeconds() - previousTime.asSeconds())); // is capped at 60
            counter = 0;
        }
        previousTime = currentTime;

        for (size_t i = 0; i < fluidSolver.numParticles; i++)
        {
            if (fluidSolver.particles[i].isFluid)
            {
                maxVelocity = max(maxVelocity, sqrt(fluidSolver.particles[i].velocity.x * fluidSolver.particles[i].velocity.x + fluidSolver.particles[i].velocity.y * fluidSolver.particles[i].velocity.y));
            }
        }

        cflNumber.setString("CFL: lambda >= " + to_string((fluidSolver.TIME_STEP * maxVelocity) / fluidSolver.H) + ", maxTimeStep: " + to_string(fluidSolver.H / maxVelocity) + ", fps: " + to_string(fps));
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
    delete[] boundingBox;

    return EXIT_SUCCESS;
}