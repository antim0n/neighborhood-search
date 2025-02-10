#include "zIndexSort.h"

vector<int> getParticleIndicesZI;
Handle* sortedIndicesZI = nullptr;
int maxValZI = 100;
int globalCounterZI = 100;


void zIndexSortConstruction(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // compute cell index with z-curve / morton code / bit interleaving
    int biggerXY = boundingBox[5];
    if (boundingBox[4] > boundingBox[5])
    {
        biggerXY = boundingBox[4];
    }
    biggerXY = pow(2, static_cast<int>(log2(biggerXY) + 1));
    int requiredSize = biggerXY * biggerXY; // z-index can be alot bigger -> always in 2^n rekursive boxes, unnessesary used memory TODO: compression
    getParticleIndicesZI.resize(requiredSize); // every cell points to the first of its particles
    fill(getParticleIndicesZI.begin(), getParticleIndicesZI.end(), 0);

    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h)); // conversion fails sometimes on edge cases like 1.000 -> 0
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));

        particles[i].k = k;
        particles[i].l = l;

        uint32_t x = k;
        uint32_t y = l;
        uint64_t zIndex = interleaveBits(x, y);
        particles[i].cellIndex = zIndex;

        // increment cellIndices
        getParticleIndicesZI.at(particles[i].cellIndex) += 1;
    }
    // accumulate cellIndices
    for (size_t i = 1; i < requiredSize; i++)
    {
        getParticleIndicesZI.at(i) += getParticleIndicesZI.at(i - 1);
    }

    // sort particles with respect to their index with the help of cellIndices
    // insertion sort
    // getParticleIndices.at(0) -= 1; why did i do that?
    for (size_t i = 1; i < numParticles; i++)
    {
        Particle current = particles[i];
        getParticleIndicesZI.at(current.cellIndex) -= 1;
        int j = i - 1;
        while (j >= 0 && current.cellIndex < particles[j].cellIndex)
        {
            particles[j + 1] = particles[j];
            j -= 1;
        }
        particles[j + 1] = current;
    }
}

void zIndexSortConstructionHandleSort(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // compute cell index with z-curve / morton code / bit interleaving
    int biggerXY = boundingBox[5];
    if (boundingBox[4] > boundingBox[5])
    {
        biggerXY = boundingBox[4];
    }
    biggerXY = pow(2, static_cast<int>(log2(biggerXY) + 1));
    int requiredSize = biggerXY * biggerXY; // z-index can be alot bigger -> always in 2^n rekursive boxes, unnessesary used memory TODO: compression
    getParticleIndicesZI.resize(requiredSize); // every cell points to the first of its particles
    fill(getParticleIndicesZI.begin(), getParticleIndicesZI.end(), 0);

    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h)); // conversion fails sometimes on edge cases like 1.000 -> 0
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));

        particles[i].k = k;
        particles[i].l = l;

        uint32_t x = k;
        uint32_t y = l;
        uint64_t zIndex = interleaveBits(x, y);
        particles[i].cellIndex = zIndex;

        // increment cellIndices
        getParticleIndicesZI.at(particles[i].cellIndex) += 1;
    }
    // accumulate cellIndices
    for (size_t i = 1; i < requiredSize; i++)
    {
        getParticleIndicesZI.at(i) += getParticleIndicesZI.at(i - 1);
    }

    // sort particles with respect to their index with the help of cellIndices
    // array with only current index and cellIndex
    if (sortedIndicesZI == nullptr)
    {
        sortedIndicesZI = new Handle[numParticles];
    }
    if (globalCounterZI == 0)
    {
        for (size_t i = 1; i < numParticles; i++)
        {
            Particle current = particles[i];
            int j = i - 1;
            while (j >= 0 && current.cellIndex < particles[j].cellIndex)
            {
                particles[j + 1] = particles[j];
                j -= 1;
            }
            particles[j + 1] = current;
        }
    }
    for (size_t i = 0; i < numParticles; i++)
    {
        sortedIndicesZI[i] = { particles[i].cellIndex, &particles[i] };
    }
    // insertion sort
    for (size_t i = 1; i < numParticles; i++)
    {
        Handle current = sortedIndicesZI[i];
        getParticleIndicesZI.at(current.cellIndex) -= 1;
        int j = i - 1;
        while (j >= 0 && current.cellIndex < sortedIndicesZI[j].cellIndex)
        {
            sortedIndicesZI[j + 1] = sortedIndicesZI[j];
            j -= 1;
        }
        sortedIndicesZI[j + 1] = current;
    }
    // initial sort
    if (globalCounterZI == maxValZI)
    {
        // copy particles twice
        Particle* sortedParticles = new Particle[numParticles];
        for (size_t i = 0; i < numParticles; i++)
        {
            sortedParticles[i] = *sortedIndicesZI[i].reference;
        }
        copy(sortedParticles, sortedParticles + numParticles, particles);
        delete[] sortedParticles;
        for (size_t i = 0; i < numParticles; i++)
        {
            sortedIndicesZI[i] = { particles[i].cellIndex, &particles[i] };
        }
    }
    globalCounterZI += 1;
    globalCounterZI %= maxValZI; // only sort particle data every nth step
}

void zIndexSortQuery(Particle* particles, int numParticles, float h)
{
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            particles[i].neighbors.clear();

            // mask to compute adjacent cell indices
            int cellIndices[][2] = {
                {-1, 1}, {0, 1}, {1, 1},
                {-1, 0}, {0, 0}, {1, 0},
                {-1, -1}, {0, -1}, {1, -1}
            };

            // check particles in all adjacent cells
            for (size_t j = 0; j < 9; j++)
            {
                int newIndexX = particles[i].k + cellIndices[j][0];
                int newIndexY = particles[i].l + cellIndices[j][1];

                // check for out of bounds cells
                if (newIndexX >= 0 && newIndexX < boundingBox[4] && newIndexY >= 0 && newIndexY < boundingBox[5])
                {
                    // cell z-index
                    uint32_t x = newIndexX;
                    uint32_t y = newIndexY;
                    uint64_t zIndex = interleaveBits(x, y);

                    // compare with current particle
                    for (size_t k = getParticleIndicesZI.at(zIndex); k < getParticleIndicesZI.at(zIndex + 1); k++)
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

void zIndexSortQueryHandleSort(Particle* particles, int numParticles, float h)
{
    for (size_t i = 0; i < numParticles; i++)
    {
        if (sortedIndicesZI[i].reference->isFluid)
        {
            sortedIndicesZI[i].reference->neighbors.clear();

            // mask to compute adjacent cell indices
            int cellIndices[][2] = {
                {-1, 1}, {0, 1}, {1, 1},
                {-1, 0}, {0, 0}, {1, 0},
                {-1, -1}, {0, -1}, {1, -1}
            };

            // check particles in all adjacent cells
            for (size_t j = 0; j < 9; j++)
            {
                int newIndexX = sortedIndicesZI[i].reference->k + cellIndices[j][0];
                int newIndexY = sortedIndicesZI[i].reference->l + cellIndices[j][1];

                // check for out of bounds cells
                if (newIndexX >= 0 && newIndexX < boundingBox[4] && newIndexY >= 0 && newIndexY < boundingBox[5])
                {
                    // cell z-index
                    uint32_t x = newIndexX;
                    uint32_t y = newIndexY;
                    uint64_t zIndex = interleaveBits(x, y);

                    // compare with current particle
                    for (size_t k = getParticleIndicesZI.at(zIndex); k < getParticleIndicesZI.at(zIndex + 1); k++)
                    {
                        Vector2f d = sortedIndicesZI[i].reference->position - sortedIndicesZI[k].reference->position;
                        float distance = sqrt(d.x * d.x + d.y * d.y);
                        if (distance < 2.0f * h)
                        {
                            sortedIndicesZI[i].reference->neighbors.push_back(sortedIndicesZI[k].reference);
                        }
                    }
                }
            }
        }
    }
}