#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>

using namespace sf;
using namespace std;

const double PI = 3.14159265358979323846;

/* holds all important information needed at a particle */
struct Particle
{
    int index;
    float density;  // in kg / m^3, 997 for water
    float pressure; // in N/m^2
    float mass;  // density * volume
    Vector2f position;
    Vector2f velocity;
    Vector2f acceleration;
    vector<Particle*> neighbors;    // pointers to all current neighbors
};

class FluidSolver // add access
{
public:
    const float H = 0.025f;
    const float REST_DENSITY = 1.2f;
    const float PRESSURE = 0.f;
    const Vector2f GRAVITY = Vector2f(0.f, -9.81f);

    const float STIFFNESS = 300.f; //
    const float VISCOSITY = 0.02f; //
    const float TIME_STEP = 0.0012f; //

    int numFluidParticles;
    int numBoundaryParticles = 380;
    int numParticles = numFluidParticles + numBoundaryParticles;

    // constructor
    FluidSolver(int size)
    {
        numFluidParticles = size;
        numParticles = numBoundaryParticles + size;
    }

    // destructor
    ~FluidSolver()
    {

    }

    // function declarations
    /* */
    void initializeFluidParticles(Particle* particles, Vector2f offset);
    /* */
    void initializeBoundaryParticles(Particle* particles);
    /* */
    float cubicSpline(Vector2f positionA, Vector2f positionB);
    /* */
    Vector2f cubicSplineDerivative(Vector2f positionA, Vector2f positionB);
    /* */
    void neighborSearchNN(Particle* particles, float support);
    /* */
    void computeDensityAndPressure(Particle* particles);
    /* */
    void updatePositions(Particle* particles);
    /* */
    void computeAccelerations(Particle* particles);
    /* */
    Vector2f nonPressureAcceleration(Particle p);
    /* */
    Vector2f pressureAcceleration(Particle p);
};
