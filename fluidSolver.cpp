#include "fluidSolver.h"
#include <iostream>


FluidSolver::FluidSolver(int size)
{
    numFluidParticles = size;
    numBoundaryParticles = sqrt(numFluidParticles) * 16;
    numParticles = numBoundaryParticles + numFluidParticles;
    particles = new Particle[numParticles];
}

FluidSolver::~FluidSolver()
{
    delete[] particles;
}

void FluidSolver::initializeFluidParticles(Vector2f offset)
{
    int size = static_cast<int>(sqrt(numFluidParticles));
    for (size_t i = 0; i < numFluidParticles; i++)
    {
        particles[i].cellIndex = -1;
        particles[i].isFluid = true;
        particles[i].density = REST_DENSITY;
        particles[i].pressure = PRESSURE;
        particles[i].mass = REST_DENSITY * H * H;
        particles[i].velocity = Vector2f(0, 0);
        particles[i].acceleration = Vector2f(0, 0);

        float temp1 = H * (i % size) + offset.x; // use a 100x100 global space
        float temp2 = H * (i / size) + offset.y;
        particles[i].index = i;
        particles[i].position = Vector2f(temp1, temp2); // distribute the particles
    }
}

void FluidSolver::initializeBoundaryParticles()
{
    for (size_t i = numFluidParticles; i < numParticles; i++)
    {
        float temp1 = 0.f;
        float temp2 = 0.f;

        int fourth = numBoundaryParticles / 4;
        if (fourth % 2)
        {
            fourth -= 1;
        }

        // basic rectangle
        if (i - numFluidParticles < fourth) // 120
        {
            // right
            temp1 = H * ((i - numFluidParticles) % 2 + fourth / 2);
            temp2 = H * ((i - numFluidParticles) / 2 + 5);
        }
        else if (i - numFluidParticles < fourth * 2)
        {
            // left
            temp1 = H * ((i - numFluidParticles - fourth) % 2 + 2);
            temp2 = H * ((i - numFluidParticles - fourth) / 2 + 5);
        }
        else if (i - numFluidParticles < fourth * 3)
        {
            // top
            temp1 = H * ((i - numFluidParticles - fourth * 2) % (fourth / 2) + 2);
            temp2 = H * ((i - numFluidParticles - fourth * 2) / (fourth / 2) + fourth / 2 + 5);
        }
        else if (i - numFluidParticles < fourth * 4)
        {
            // bottom
            temp1 = H * ((i - numFluidParticles - fourth * 3) % (fourth / 2) + 2);
            temp2 = H * ((i - numFluidParticles - fourth * 3) / (fourth / 2) + 3);
        }
        else
        {
            // bottom
            temp1 = H * ((i - numFluidParticles - fourth * 4) % (fourth / 2) + 2 + fourth / 2);
            temp2 = H * ((i - numFluidParticles - fourth * 4) / (fourth / 2) + 3);
        }

        particles[i].pressure = PRESSURE;
        particles[i].velocity = Vector2f(0, 0);
        particles[i].acceleration = Vector2f(0, 0);

        particles[i].index = i;
        particles[i].cellIndex = -1;
        particles[i].isFluid = false;
        particles[i].position = Vector2f(temp1, temp2); // distribute the particles
        particles[i].density = REST_DENSITY;
        particles[i].mass = REST_DENSITY * H * H;
    }
}

float FluidSolver::cubicSpline(Vector2f positionA, Vector2f positionB) const // just slightly too much * 0.999138886f for alpha
{
    Vector2f temp = positionA - positionB;
    float distance = sqrt(temp.x * temp.x + temp.y * temp.y);
    float d = distance / H;

    float alpha = 5.f / (14.f * PI * H * H); // 2D

    float t1 = max(1.f - d, 0.f);
    float t2 = max(2.f - d, 0.f);

    float w = alpha * (t2 * t2 * t2 - 4.f * t1 * t1 * t1);

    return w;
}

Vector2f FluidSolver::cubicSplineDerivative(Vector2f positionA, Vector2f positionB)
{
    Vector2f temp = positionA - positionB;
    float distance = sqrt(temp.x * temp.x + temp.y * temp.y);
    if (distance == 0.f)
    {
        return Vector2f(0.f, 0.f);
    }
    float d = distance / H;

    float alpha = 5.f / (14.f * PI * H * H); // 2D

    float t1 = max(1 - d, 0.f);
    float t2 = max(2 - d, 0.f);

    float v = alpha * (-3.f * t2 * t2 + 12.f * t1 * t1);
    Vector2f w1 = Vector2f((temp.x / (d * H * H)) * v, (temp.y / (d * H * H)) * v); // more accurate kernel wiht a scalar?

    return w1;
}

void FluidSolver::neighborSearchNN(float support)
{
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            particles[i].neighbors.clear();
            for (size_t j = 0; j < numParticles; j++)
            {
                Vector2f d = particles[i].position - particles[j].position;
                float distance = sqrt(d.x * d.x + d.y * d.y);
                if (distance < support * H)
                {
                    particles[i].neighbors.push_back(&particles[j]);
                }
            }
        }
    }
}

void FluidSolver::computeDensityAndPressure()
{
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            float temp = 0;
            if (particles[i].neighbors.size() == 1)
            {
                temp = REST_DENSITY;
            }
            else
            {
                // sum over all neighbors
                for (size_t j = 0; j < particles[i].neighbors.size(); j++)
                {
                    temp += particles[i].neighbors[j]->mass * cubicSpline(particles[i].position, particles[i].neighbors[j]->position);
                }
            }
            particles[i].density = temp;
            particles[i].pressure = STIFFNESS * (max((temp / REST_DENSITY) - 1.f, 0.f)); // problem with dividing by 0
        }
    }
}

void FluidSolver::updatePositions()
{
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            particles[i].velocity += TIME_STEP * particles[i].acceleration;
            particles[i].position += TIME_STEP * particles[i].velocity; // updated velocity (semi-implicit euler)
        }
    }
}

Vector2f FluidSolver::nonPressureAcceleration(Particle p)
{
    Vector2f SPH = Vector2f(0.f, 0.f);
    for (size_t i = 0; i < p.neighbors.size(); i++)
    {
        Vector2f velD = p.velocity - p.neighbors[i]->velocity;
        Vector2f posD = p.position - p.neighbors[i]->position;
        float val = (velD.x * posD.x + velD.y * posD.y) / (posD.x * posD.x + posD.y * posD.y + 0.01f * H * H);
        Vector2f kernel = cubicSplineDerivative(p.position, p.neighbors[i]->position);
        SPH += (p.neighbors[i]->mass / p.neighbors[i]->density) * val * kernel;
    }
    return 2.f * VISCOSITY * SPH + GRAVITY;
}

Vector2f FluidSolver::pressureAcceleration(Particle p)
{
    Vector2f SPH = Vector2f(0.f, 0.f);
    for (size_t i = 0; i < p.neighbors.size(); i++)
    {
        float val = 0;
        if (!p.neighbors[i]->isFluid) // boundary handling (mirroring)
        {
            val = p.pressure / (p.density * p.density) + p.pressure / (p.density * p.density);
        }
        else
        {
            val = p.pressure / (p.density * p.density) + p.neighbors[i]->pressure / (p.neighbors[i]->density * p.neighbors[i]->density);
        }
        Vector2f kernel = cubicSplineDerivative(p.position, p.neighbors[i]->position);
        SPH += p.neighbors[i]->mass * val * kernel;
    }
    return -SPH;
}

void FluidSolver::computeAccelerations()
{
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            Vector2f aNonP = nonPressureAcceleration(particles[i]);
            Vector2f aP = pressureAcceleration(particles[i]);

            particles[i].acceleration = aNonP + aP;
        }
    }
}