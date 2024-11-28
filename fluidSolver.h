#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"

using namespace std;
using namespace sf;

const double PI = 3.14159265358979323846;

class FluidSolver // add access
{
public:
    const float STIFFNESS = 14400.f; //
    const float VISCOSITY = 7.95f; //
    const float TIME_STEP = 0.0081f; //

    const float H = 1.2f;
    const float REST_DENSITY = 1.2f;
    const float PRESSURE = 0.f;
    const Vector2f GRAVITY = Vector2f(0.f, -9.81f);

    int numFluidParticles;
    int numBoundaryParticles = 400;
    int numParticles;

    Particle* particles; // maybe change to vector

    // constructor
    FluidSolver(int size);
    // destructor
    ~FluidSolver();

    /* */
    Vector2f particleToPixelCoord(Vector2f particlePos, int windowWidth, int windowHeight);
    /* */
    Vector2f pixelToParticleCoord(Vector2f pixelPos, int windowWidth, int windowHeight);
    /* */
    void initializeFluidParticles(Vector2f offset);
    /* */
    void initializeBoundaryParticles();
    /* */
    void neighborSearchNN(float support);
    /* */
    void computeDensityAndPressure();
    /* */
    void updatePositions();
    /* */
    void computeAccelerations();
    /* */
    float cubicSpline(Vector2f positionA, Vector2f positionB) const;
    /* */
    Vector2f cubicSplineDerivative(Vector2f positionA, Vector2f positionB);
    /* */
    Vector2f nonPressureAcceleration(Particle p);
    /* */
    Vector2f pressureAcceleration(Particle p);
};
