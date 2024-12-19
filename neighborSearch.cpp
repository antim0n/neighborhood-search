#include "neighborSearch.h"
#include <iostream>
#include <chrono>

using namespace std::chrono;

vector<vector<Particle*>> cells = vector<vector<Particle*>>();
float* boundingBox = new float[6];
vector<int> getParticleIndices = vector<int>();

void boundingBoxConstruction(Particle* particles, int numParticles, float h)
{
    // xmin, xmax, ymin, ymax, cellsx, cellsy
    boundingBox[0] = particles[0].position.x; // min
    boundingBox[1] = particles[0].position.x; // max
    boundingBox[2] = particles[0].position.y; // min
    boundingBox[3] = particles[0].position.y; // max
    // compute bounding box
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].position.x < boundingBox[0])
        {
            boundingBox[0] = particles[i].position.x;
        }
        else if (particles[i].position.x > boundingBox[1])
        {
            boundingBox[1] = particles[i].position.x;
        }
        if (particles[i].position.y < boundingBox[2])
        {
            boundingBox[2] = particles[i].position.y;
        }
        else if (particles[i].position.y > boundingBox[3])
        {
            boundingBox[3] = particles[i].position.y;
        }
    }
    boundingBox[0] -= 0.4 * h; // make grid slightly bigger for edge cases
    boundingBox[1] += 0.4 * h;
    boundingBox[2] -= 0.4 * h;
    boundingBox[3] += 0.4 * h;

    boundingBox[4] = ceil((boundingBox[1] - boundingBox[0]) / (2.f * h));
    boundingBox[5] = ceil((boundingBox[3] - boundingBox[2]) / (2.f * h));
}

void gridConstruction(Particle* particles, int numParticles, float h)
{
    // uniform grid
    // particles are stored in a cell, without sorting
    boundingBoxConstruction(particles, numParticles, h);

    // compute cell index with (k, l, m)
    for (size_t i = 0; i < boundingBox[4] * boundingBox[5] + 1; i++)
    {
        cells.push_back(vector<Particle*>());
    }
    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h));
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));
        particles[i].cellIndex = k + l * boundingBox[4];
        cells.at(particles[i].cellIndex).push_back(&particles[i]);
    }
}

void gridQuery(Particle* particles, int numParticles, float h) // wrong
{
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            particles[i].neighbors.clear();
            // compute cell indices
            int cellIndices[] = { particles[i].cellIndex,
                particles[i].cellIndex + 1,
                particles[i].cellIndex - 1,
                particles[i].cellIndex + boundingBox[4],
                particles[i].cellIndex - boundingBox[4],
                particles[i].cellIndex + boundingBox[4] + 1,
                particles[i].cellIndex + boundingBox[4] - 1,
                particles[i].cellIndex - boundingBox[4] + 1,
                particles[i].cellIndex - boundingBox[4] - 1
            };

            // check particles in all adjacent cells
            for (size_t j = 0; j < 9; j++)
            {
                if (cellIndices[j] >= 0 && cellIndices[j] < boundingBox[4] * boundingBox[5]) // valid cell index
                {
                    for (size_t k = 0; k < cells.at(cellIndices[j]).size(); k++)
                    {
                        Vector2f d = particles[i].position - cells.at(cellIndices[j]).at(k)->position;
                        float distance = sqrt(d.x * d.x + d.y * d.y);
                        if (distance < 2.0f * h)
                        {
                            particles[i].neighbors.push_back(cells.at(cellIndices[j]).at(k));
                        }
                    }
                }
            }
        }
    }
}

// index sort takes particles as input, sorts them and returns an integer array as output
void indexSortConstruction(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // compute cell index with (k, l, m)
    getParticleIndices.clear();
    for (size_t i = 0; i < boundingBox[4] * boundingBox[5] + 1; i++)
    {
        getParticleIndices.push_back(0);
    }
    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h)); // conversion fails sometimes on edge cases like 1.000 -> 0
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));
        particles[i].cellIndex = k + l * boundingBox[4]; // TODO: with morton code (z-curve, bit-interleaving)
        // increment cellIndices
        getParticleIndices.at(particles[i].cellIndex) += 1;
    }
    // accumulate cellIndices
    for (size_t i = 1; i < boundingBox[4] * boundingBox[5] + 1; i++)
    {
        getParticleIndices.at(i) += getParticleIndices.at(i - 1);
    }
    // sort particles with respect to their index with the help of cellIndices
    // insertion sort
    getParticleIndices.at(0) -= 1;
    for (size_t i = 1; i < numParticles; i++)
    {
        Particle current = particles[i];
        getParticleIndices.at(current.cellIndex) -= 1;
        int j = i - 1;
        while (j >= 0 && current.cellIndex < particles[j].cellIndex)
        {
            particles[j + 1] = particles[j];
            j -= 1;
        }
        particles[j + 1] = current;
    }
}

void indexSortQuery(Particle* particles, int numParticles, float h) // no boundary particles included, too many adjacent cells at times, maybe check with klm
{
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            particles[i].neighbors.clear();
            // compute cell indices
            int cellIndices[] = { particles[i].cellIndex,
                particles[i].cellIndex + 1,
                particles[i].cellIndex - 1,
                particles[i].cellIndex + boundingBox[4],
                particles[i].cellIndex - boundingBox[4],
                particles[i].cellIndex + boundingBox[4] + 1,
                particles[i].cellIndex + boundingBox[4] - 1,
                particles[i].cellIndex - boundingBox[4] + 1,
                particles[i].cellIndex - boundingBox[4] - 1
            };

            // check particles in all adjacent cells
            for (size_t j = 0; j < 9; j++)
            {
                if (cellIndices[j] >= 0 && cellIndices[j] < boundingBox[4] * boundingBox[5]) // valid cell index
                {
                    for (size_t k = getParticleIndices.at(cellIndices[j]); k < getParticleIndices.at(cellIndices[j] + 1); k++)
                    {
                        Vector2f d = particles[i].position - particles[k].position;
                        float distance = sqrt(d.x * d.x + d.y * d.y);
                        if (distance < 2.0f * h)
                        {
                            particles[i].neighbors.push_back(&particles[k]);
                        }
                    }
                }
            }
        }
    }
}

void zIndexSortConstruction()
{

}

void zIndexSortQuery()
{

}

void spatialHashingConstruction()
{
    // maps grid cell to a hash cell
    // compute cell index c or cell identifier (k, l, m) for particles
    // compute hash function i = h(c) or i = h(k, l, m)
    // store particles in array (hash table) at index i (maybe vector of vectors?)
}
void spatialHashingQuery()
{
    // with hash function h(c)
}

void compactHashingConstruction()
{

}

void compactHashingQuery()
{

}