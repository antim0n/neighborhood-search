#include "indexSort.h"


vector<int> getParticleIndicesI;

// index sort takes particles as input, sorts them and returns an integer array as output
void indexSortConstruction(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // compute cell index with (k, l, m)
    getParticleIndicesI.resize(boundingBox[4] * boundingBox[5] + 1);
    fill(getParticleIndicesI.begin(), getParticleIndicesI.end(), 0);

    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h)); // conversion fails sometimes on edge cases like 1.000 -> 0
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));
        particles[i].cellIndex = k + l * boundingBox[4];
        // increment cellIndices
        getParticleIndicesI.at(particles[i].cellIndex) += 1;
    }
    // accumulate cellIndices
    for (size_t i = 1; i < boundingBox[4] * boundingBox[5] + 1; i++)
    {
        getParticleIndicesI.at(i) += getParticleIndicesI.at(i - 1);
    }
    // sort particles with respect to their index with the help of cellIndices
    // insertion sort
    // getParticleIndices.at(0) -= 1;
    for (size_t i = 1; i < numParticles; i++)
    {
        Particle current = particles[i];
        getParticleIndicesI.at(current.cellIndex) -= 1;
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
                if (cellIndices[j] >= 0 && cellIndices[j] < boundingBox[4] * boundingBox[5]) // valid cell index TODO: change to k,l
                {
                    for (size_t k = getParticleIndicesI.at(cellIndices[j]); k < getParticleIndicesI.at(cellIndices[j] + 1); k++)
                    {
                        Vector2f d = particles[i].position - particles[k].position;
                        float distance = d.x * d.x + d.y * d.y;
                        if (distance < (2.0f * h) * (2.0f * h))
                        {
                            particles[i].neighbors.push_back(&particles[k]);
                        }
                    }
                }
            }
        }
    }
}