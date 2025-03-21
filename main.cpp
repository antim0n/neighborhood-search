#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include <iostream>
#include <chrono>
#include <numeric>
#include "button.h"
#include "fluidSolver.h"
#include "basicGrid.h"
#include "indexSort.h"
#include "zIndexSort.h"
#include "spatialHashing.h"
#include "compactHashing.h"
#include "cellLinkedList.h"

using namespace std;
using namespace std::chrono;
using namespace sf;


static Vector2f particleToPixelCoord(Vector2f particlePos, int windowWidth, int windowHeight, float scaling)
{
    return Vector2f((particlePos.x) * (static_cast<float>(windowWidth) / scaling),
        static_cast<float>(windowHeight) - particlePos.y * (static_cast<float>(windowWidth) / scaling));
}

static Vector2f pixelToParticleCoord(Vector2f pixelPos, int windowWidth, int windowHeight, float scaling)
{
    return Vector2f(static_cast<float>(pixelPos.x) / (static_cast<float>(windowWidth) / scaling),
        static_cast<float>(windowHeight - pixelPos.y) / (static_cast<float>(windowWidth) / scaling));
}

int main()
{
    FluidSolver fluidSolver(1000);

    /* initialize all particles */
    fluidSolver.initializeFluidParticles(Vector2f(sqrt(fluidSolver.numFluidParticles) / 2 * fluidSolver.H, 5.f * fluidSolver.H));
    fluidSolver.initializeBoundaryParticles();
    int scaling = sqrt(fluidSolver.numFluidParticles) * 3; // bigger -> smaller picture

    /*fluidSolver.LoadMaze();
    int scaling = sqrt(fluidSolver.numBoundaryParticles) * 12;*/

    const int WINDOW_WIDTH = 900;
    const int WINDOW_HEIGHT = 900;

    // states
    bool stopSimulation = false;
    bool showNeighbors = false;
    bool showGrid = false;
    bool gatherRuntimes = false;

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

    // performance average
    vector<chrono::milliseconds> runtimeConstruction;
    vector<chrono::milliseconds> runtimeQuery;
    vector<chrono::milliseconds> runtimeTotal;

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
    Button b1(Vector2i(50, 20), Vector2i(10, 50), Color::Green, Text("grid", font, 15));
    Button b2(Vector2i(60, 20), Vector2i(70, 50), Color::Green, Text("runtime", font, 15));
    Button b3(Vector2i(70, 20), Vector2i(140, 50), Color::Green, Text("maze", font, 15));

    /* allocate memory for the particle shapes */
    int numVertices = 12;
    VertexArray* drawingCircles = new VertexArray[fluidSolver.numParticles];
    float toRadian = PI / 180.f;
    for (size_t i = 0; i < fluidSolver.numParticles; i++)
    {
        VertexArray newCircle = VertexArray(TriangleFan, numVertices);
        newCircle[0].position = Vector2f(0.f, 0.f);
        newCircle[1].position = Vector2f(1.f * cos(0.f * toRadian), 1.f * sin(0.f * toRadian));
        newCircle[2].position = Vector2f(1.f * cos(36.f * toRadian), 1.f * sin(36.f * toRadian));
        newCircle[3].position = Vector2f(1.f * cos(72.f * toRadian), 1.f * sin(72.f * toRadian));
        newCircle[4].position = Vector2f(1.f * cos(108.f * toRadian), 1.f * sin(108.f * toRadian));
        newCircle[5].position = Vector2f(1.f * cos(144.f * toRadian), 1.f * sin(144.f * toRadian));
        newCircle[6].position = Vector2f(1.f * cos(180.f * toRadian), 1.f * sin(180.f * toRadian));
        newCircle[7].position = Vector2f(1.f * cos(216.f * toRadian), 1.f * sin(216.f * toRadian));
        newCircle[8].position = Vector2f(1.f * cos(252.f * toRadian), 1.f * sin(252.f * toRadian));
        newCircle[9].position = Vector2f(1.f * cos(288.f * toRadian), 1.f * sin(288.f * toRadian));
        newCircle[10].position = Vector2f(1.f * cos(324.f * toRadian), 1.f * sin(324.f * toRadian));
        newCircle[11].position = Vector2f(1.f * cos(360.f * toRadian), 1.f * sin(360.f * toRadian));
        drawingCircles[i] = newCircle;
    }
    if (!fluidSolver.particles)
    {
        cout << "Memory allocation failed.\n";
    }

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
                if (event.mouseButton.button == sf::Mouse::Right) // TODO: fix, doesnt work because sorted
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
                        showGrid = !showGrid;
                    }
                    if (b2.border.contains(Vector2i(Mouse::getPosition(window))))
                    {
                        if (gatherRuntimes)
                        {
                            cout << "End collecting ------------------------" << endl;

                            auto duration = accumulate(runtimeConstruction.begin(), runtimeConstruction.end(), decltype(runtimeConstruction)::value_type(0));
                            duration = duration / runtimeConstruction.size();
                            cout << "Construction: " << duration.count() << " " << chrono::duration<double>(duration).count() << endl;

                            auto duration2 = accumulate(runtimeQuery.begin(), runtimeQuery.end(), decltype(runtimeQuery)::value_type(0));
                            duration2 = duration2 / runtimeQuery.size();
                            cout << "Query: " << duration2.count() << " " << chrono::duration<double>(duration2).count() << endl;

                            auto duration3 = accumulate(runtimeTotal.begin(), runtimeTotal.end(), decltype(runtimeTotal)::value_type(0));
                            duration3 = duration3 / runtimeTotal.size();
                            cout << "Total: " << duration3.count() << " " << chrono::duration<double>(duration3).count() << endl;
                            cout << runtimeTotal.size() << endl;
                            cout << endl;

                            gatherRuntimes = false;
                        }
                        else
                        {
                            cout << "Start collecting ------------------------" << endl;
                            runtimeConstruction.clear();
                            runtimeQuery.clear();
                            runtimeTotal.clear();
                            gatherRuntimes = true;
                        } 
                    }
                    if (b3.border.contains(Vector2i(Mouse::getPosition(window))))
                    {
                        // fluidSolver.LoadMaze(); // e.g. circle shapes would also have to be reallocated
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

            // boundingBoxConstruction(); TODO
            // initialize(); precompute some global variable depending on the algorithm

            fluidSolver.neighborSearchNN(2);
            
            // gridConstruction(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // gridConstructionImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // gridConstructionImprovedParallel(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            
            // indexSortConstructionCountingSort(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // indexSortConstructionCompareSorting(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // indexSortConstructionInsertionSortImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // indexSortConstructionCountingSortImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // indexSortConstructionCountingSortImprovedParallel(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);

            // zIndexSortConstruction(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortConstructionImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortConstructionHandleSort(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortConstructionHandleSortImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortConstructionHandleSortImprovedMap(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortConstructionHandleSortImprovedParallel(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            
            // spatialHashingConstruction(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // spatialHashingConstructionImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // spatialHashingConstructionImprovedParallel(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);

            // compactHashingConstruction(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // compactHashingConstructionImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // compactHashingConstructionHashCollisionFlagImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // compactHashingConstructionZSortedImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // compactHashingConstructionHandleSort(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // compactHashingConstructionHandleSortImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // compactHashingConstructionHandleSortImprovedParallel(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);

            // cellLinkedListConstruction(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);

            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<chrono::milliseconds>(stop - start);
            cout << "Construction: " << duration.count() << " " << chrono::duration<double>(duration).count() << endl;

            auto start2 = high_resolution_clock::now();

            // gridQuery(fluidSolver.particles, fluidSolver.numFluidParticles, fluidSolver.H);
            // gridQueryOverCells(fluidSolver.particles, fluidSolver.H);
            // gridQueryImproved(fluidSolver.particles, fluidSolver.numFluidParticles, fluidSolver.H);
            // gridQueryImprovedParallel(fluidSolver.particles, fluidSolver.numFluidParticles, fluidSolver.H);

            // indexSortQueryOverCells(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // indexSortQueryOverCellsImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // indexSortQuery(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // indexSortQueryKLImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // indexSortQueryImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // indexSortQueryImprovedParallel(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);

            // zIndexSortQuery(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortQueryImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortQueryHandleSort(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortQueryHandleSortImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortQueryHandleSortOverCellsImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortQueryHandleSortImprovedMap(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // zIndexSortQueryHandleSortImprovedParallel(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            
            // spatialHashingQuery(fluidSolver.particles, fluidSolver.numFluidParticles, fluidSolver.H);
            // spatialHashingQueryImproved(fluidSolver.particles, fluidSolver.numFluidParticles, fluidSolver.H);
            // spatialHashingQueryImprovedParallel(fluidSolver.particles, fluidSolver.numFluidParticles, fluidSolver.H);

            // compactHashingQuery(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // compactHashingQueryImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // compactHashingQueryHashCollisionFlagImproved(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);
            // compactHashingQueryImprovedParallel(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);

            // cellLinkedListQuery(fluidSolver.particles, fluidSolver.numParticles, fluidSolver.H);

            /*for (size_t i = 0; i < fluidSolver.numFluidParticles; i++)
            {
                for (size_t j = 0; j < fluidSolver.particles[i].neighbors.size(); j++)
                {
                    cout << fluidSolver.particles[fluidSolver.particles[i].neighbors.at(j)].index << ", ";
                }
                cout << endl;
            }*/

            auto stop2 = high_resolution_clock::now();
            auto duration2 = duration_cast<chrono::milliseconds>(stop2 - start2);
            cout << "Query: " << duration2.count() << " " << chrono::duration<double>(duration2).count() << endl;

            auto start3 = high_resolution_clock::now();

            fluidSolver.computeDensityAndPressure();
            fluidSolver.computeAccelerations();
            // fluidSolver.computeDensityAndPressureCLL();
            // fluidSolver.computeAccelerationsCLL();
            fluidSolver.updatePositions();

            // stopSimulation = true;

            auto stop3 = high_resolution_clock::now();
            // auto duration3 = duration_cast<chrono::milliseconds>(stop3 - start3);
            // cout << "Physics: " << duration3.count() << " " << chrono::duration<double>(duration3).count() << endl;
            auto duration4 = duration_cast<chrono::milliseconds>(stop3 - start);
            cout << "Total: " << duration4.count() << " " << chrono::duration<double>(duration4).count() << endl;
            cout << endl;

            if (gatherRuntimes)
            {
                runtimeConstruction.push_back(duration);
                runtimeQuery.push_back(duration2);
                runtimeTotal.push_back(duration4);
            }
        }

        /* Draw */

        window.setView(view);
        window.clear(); // don't draw on top of the previous frame

        for (size_t i = 0; i < fluidSolver.numParticles; i++)
        {
            Transform particleTransform;

            float radius = (fluidSolver.H / 2.f) * (static_cast<float>(WINDOW_HEIGHT) / scaling);
            if (radius < 1.f)
            {
                radius = 1.f;
            }
            particleTransform.translate(particleToPixelCoord(fluidSolver.particles[i].position, WINDOW_WIDTH, WINDOW_HEIGHT, scaling)); // the shapes to be drawn have to be updated independently, scale
            particleTransform.scale({ radius, radius }); // h is defined as the "diameter", shouldnt be smaller than 1
            if (fluidSolver.particles[i].isFluid)
            {
                for (size_t j = 0; j < numVertices; j++)
                {
                    drawingCircles[i][j].color = Color::Blue;
                }
            }
            else
            {
                for (size_t j = 0; j < numVertices; j++)
                {
                    drawingCircles[i][j].color = Color::White;
                }
            }
            window.draw(drawingCircles[i], particleTransform);
        }

        // grid
        if (showGrid)
        {
            for (size_t i = 0; i < boundingBox[4] + 1; i++) // vertical
            {
                Vertex line[] = {
                    Vertex(particleToPixelCoord(Vector2f(boundingBox[0] + 2.f * fluidSolver.H * i, boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Yellow),
                    Vertex(particleToPixelCoord(Vector2f(boundingBox[0] + 2.f * fluidSolver.H * i, boundingBox[2] + 2.f * fluidSolver.H * (boundingBox[5])), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Yellow) };
                window.draw(line, 2, Lines);
            }
            for (size_t i = 0; i < boundingBox[5] + 1; i++) // horizontal
            {
                Vertex line[] = {
                    Vertex(particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[2] + 2.f * fluidSolver.H * i), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Yellow),
                    Vertex(particleToPixelCoord(Vector2f(boundingBox[0] + 2.f * fluidSolver.H * boundingBox[4], boundingBox[2] + 2.f * fluidSolver.H * i), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Yellow) };
                window.draw(line, 2, Lines);
            }

            // bounding box
            Vertex line1[] = {
                Vertex(particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Red),
                Vertex(particleToPixelCoord(Vector2f(boundingBox[1], boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Red) };
            Vertex line2[] = {
                Vertex(particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Red),
                Vertex(particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[3]), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Red) };
            Vertex line3[] = {
                Vertex(particleToPixelCoord(Vector2f(boundingBox[0], boundingBox[3]), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Red),
                Vertex(particleToPixelCoord(Vector2f(boundingBox[1], boundingBox[3]), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Red) };
            Vertex line4[] = {
                Vertex(particleToPixelCoord(Vector2f(boundingBox[1], boundingBox[2]), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Red),
                Vertex(particleToPixelCoord(Vector2f(boundingBox[1], boundingBox[3]), WINDOW_WIDTH, WINDOW_HEIGHT, scaling), Color::Red) };
            window.draw(line1, 2, Lines);
            window.draw(line2, 2, Lines);
            window.draw(line3, 2, Lines);
            window.draw(line4, 2, Lines);
        }

        /* text */
        if (showNeighbors)
        {
            for (size_t i = 0; i < fluidSolver.numParticles; i++)
            {
                if (fluidSolver.particles[i].isFluid)
                {
                    particleLables[i].setString(to_string(fluidSolver.particles[i].neighbors.size()));
                    Vector2f pixelCoord = particleToPixelCoord(fluidSolver.particles[i].position, WINDOW_WIDTH, WINDOW_HEIGHT, scaling);
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
         window.draw(b2.shape);
         window.draw(b2.name);
        window.draw(b3.shape);
        window.draw(b3.name);

        /* Display */
        window.display();
    }

    /* deallocate memory */
    delete[] drawingCircles;
    delete[] particleLables;
    delete[] boundingBox;
    delete[] sortedIndicesZI;
    delete[] sortedIndicesCH;
    delete[] sortedParticlesI;
    delete[] hashTableSH;
    delete[] hashTableCH;
    delete[] numNeighbors;
    delete[] compactCellIndex;
    delete[] compactReference;
    delete[] sortedParticles;

    return EXIT_SUCCESS;
}