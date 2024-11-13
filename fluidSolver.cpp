#include "fluidSolver.h"
#include <iostream>


void FluidSolver::initializeFluidParticles(Particle* particles, Vector2f offset)
{
    for (size_t i = 0; i < numFluidParticles; i++)
    {
        particles[i].density = REST_DENSITY;
        particles[i].pressure = PRESSURE;
        particles[i].mass = REST_DENSITY * H * H;
        particles[i].velocity = Vector2f(0, 0);
        particles[i].acceleration = Vector2f(0, 0);

        float temp1 = H * (i % 25 + offset.x) - 1.f;
        float temp2 = H * (i / 25 + offset.y) - 1.f;
        particles[i].index = i;
        particles[i].position = Vector2f(temp1, temp2); // distribute the particles
    }
}

void FluidSolver::initializeBoundaryParticles(Particle* particles)
{
    for (size_t i = numFluidParticles; i < numParticles; i++)
    {
        float temp1 = 0.f;
        float temp2 = 0.f;

        // basic rectangle
        if (i - numFluidParticles < 120)
        {
            // right
            temp1 = H * ((i - numFluidParticles) % 2 + 70) - 1.f;
            temp2 = (H * ((i - numFluidParticles) / 2 + 5) - 1.f);
        }
        else if (i - numFluidParticles < 240)
        {
            // left
            temp1 = H * ((i - numFluidParticles - 120) % 2 + 2) - 1.f;
            temp2 = (H * ((i - numFluidParticles - 120) / 2 + 5) - 1.f);
        }
        else
        {
            // bottom
            temp1 = H * ((i - numFluidParticles - 240) % 70 + 2) - 1.f;
            temp2 = (H * ((i - numFluidParticles - 240) / 70 + 3) - 1.f);
        }

        particles[i].pressure = PRESSURE;
        particles[i].velocity = Vector2f(0, 0);
        particles[i].acceleration = Vector2f(0, 0);

        particles[i].index = i;
        particles[i].position = Vector2f(temp1, temp2); // distribute the particles
        particles[i].density = REST_DENSITY;
        particles[i].mass = REST_DENSITY * H * H;
    }
}

float FluidSolver::cubicSpline(Vector2f positionA, Vector2f positionB) // just slightly too much * 0.999138886f for alpha
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

void FluidSolver::neighborSearchNN(Particle* particles, float support)
{
    for (size_t i = 0; i < numFluidParticles; i++)
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

void FluidSolver::computeDensityAndPressure(Particle* particles)
{
    for (size_t i = 0; i < numFluidParticles; i++)
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

void FluidSolver::updatePositions(Particle* particles)
{
    for (size_t i = 0; i < numFluidParticles; i++)
    {
        particles[i].velocity += TIME_STEP * particles[i].acceleration;
        particles[i].position += TIME_STEP * particles[i].velocity; // updated velocity (semi-implicit euler)
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
        if (p.neighbors[i]->index >= numFluidParticles) // boundary handling (mirroring)
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

void FluidSolver::computeAccelerations(Particle* particles)
{
    for (size_t i = 0; i < numFluidParticles; i++)
    {
        Vector2f aNonP = nonPressureAcceleration(particles[i]);
        Vector2f aP = pressureAcceleration(particles[i]);

        particles[i].acceleration = aNonP + aP;
    }
}