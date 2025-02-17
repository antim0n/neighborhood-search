#include <iostream>
#include "zIndexSort.h"

vector<int> getParticleIndicesZI;
Handle* sortedIndicesZI = nullptr;
int maxValZI = 4;
int globalCounterZI = 4;
int requiredSize = 0;

void zIndexSortConstruction(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // compute cell index with z-curve / morton code / bit interleaving
    // next biggest of 2^n with respect to the longer side
    int biggerXY = boundingBox[5];
    if (boundingBox[4] > boundingBox[5])
    {
        biggerXY = boundingBox[4];
    }
    biggerXY = pow(2, static_cast<int>(log2(biggerXY) + 1));
    int requiredSize = biggerXY * biggerXY + 1; // z-index can be alot bigger -> always in 2^n rekursive boxes, unnessesary used memory TODO: compression
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

void zIndexSortConstructionImproved(Particle* particles, int numParticles, float h) // TODO building getParticleIndices may be improvable with a differnt query (save the accumulation)
{
    int oldBoundingBoxSize = boundingBox[6];
    boundingBoxConstruction(particles, numParticles, h);

    // prevent unnessecary calculations
    if (oldBoundingBoxSize != boundingBox[6])
    {
        int biggerXY = max(boundingBox[4], boundingBox[5]);
        biggerXY = pow(2, static_cast<int>(log2(biggerXY) + 1));
        requiredSize = biggerXY * biggerXY + 1;
        getParticleIndicesZI.resize(requiredSize);
    }
    memset(&getParticleIndicesZI[0], 0, requiredSize * sizeof(getParticleIndicesZI[0])); // faster than fill()


    // precompute
    float invCellSize = 1.0f / (2.f * h);
    int k;
    int l;
    for (size_t i = 0; i < numParticles; i++)
    {
        k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);

        // reduce interleave computation (not significant)
        if (k != particles[i].k || l != particles[i].l)
        {
            particles[i].k = k;
            particles[i].l = l;
            particles[i].cellIndex = interleaveBits(k, l);
        }
        getParticleIndicesZI[particles[i].cellIndex] += 1; // remove at()
    }
    for (size_t i = 1; i < requiredSize; i++)
    {
        getParticleIndicesZI[i] += getParticleIndicesZI[i - 1];
    }

    for (size_t i = 1; i < numParticles; i++)
    {
        Particle current = particles[i];
        getParticleIndicesZI[current.cellIndex] -= 1;
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

    // compute cell index with z-curve / morton code / bit interleaving TODO compressed neighbors
    int biggerXY = boundingBox[5];
    if (boundingBox[4] > boundingBox[5])
    {
        biggerXY = boundingBox[4];
    }
    biggerXY = pow(2, static_cast<int>(log2(biggerXY) + 1));
    int requiredSize = biggerXY * biggerXY + 1; // z-index can be alot bigger -> always in 2^n rekursive boxes, unnessesary used memory TODO: compression
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
    // sort the actual particles array every (100)th step
    if (globalCounterZI == 100)
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
    // build handle structure
    for (size_t i = 0; i < numParticles; i++)
    {
        sortedIndicesZI[i] = { particles[i].cellIndex, &particles[i] };
    }
    // insertion sort the handle structure
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
    globalCounterZI %= maxValZI; // only sort particle data every nth step
    globalCounterZI += 1;
}

void zIndexSortConstructionHandleSortImproved(Particle* particles, int numParticles, float h)
{
    int oldBoundingBoxSize = boundingBox[6];
    boundingBoxConstruction(particles, numParticles, h);

    // prevent unnessecary calculations
    if (oldBoundingBoxSize != boundingBox[6])
    {
        int biggerXY = max(boundingBox[4], boundingBox[5]);
        biggerXY = pow(2, static_cast<int>(log2(biggerXY) + 1));
        requiredSize = biggerXY * biggerXY + 1;
        getParticleIndicesZI.resize(requiredSize);
    }
    memset(&getParticleIndicesZI[0], 0, requiredSize * sizeof(getParticleIndicesZI[0])); // is faster than fill()


    // precompute and preallocate
    float invCellSize = 1.0f / (2.f * h);
    int k;
    int l;
    for (size_t i = 0; i < numParticles; i++)
    {
        k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);

        particles[i].k = k;
        particles[i].l = l;

        particles[i].cellIndex = interleaveBits(k, l); // reduce variables
        getParticleIndicesZI[particles[i].cellIndex] += 1; // at() to []
    }
    for (size_t i = 1; i < requiredSize; i++)
    {
        getParticleIndicesZI[i] += getParticleIndicesZI[i - 1];
    }

    if (sortedIndicesZI == nullptr)
    {
        sortedIndicesZI = new Handle[numParticles];
    }
    if (globalCounterZI == 0) // TODO maybe also copy twice
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
    for (size_t i = 1; i < numParticles; i++)
    {
        Handle current = sortedIndicesZI[i];
        getParticleIndicesZI[current.cellIndex] -= 1;
        int j = i - 1;
        while (j >= 0 && current.cellIndex < sortedIndicesZI[j].cellIndex)
        {
            sortedIndicesZI[j + 1] = sortedIndicesZI[j];
            j -= 1;
        }
        sortedIndicesZI[j + 1] = current;
    }
    // initial sort (executed only once)
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
    if (globalCounterZI >= maxValZI)
    {
        globalCounterZI = 0;
    }
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

void zIndexSortQueryImproved(Particle* particles, int numParticles, float h)
{
    // precompute and allocate
    int h2 = (2.0f * h) * (2.0f * h);
    int newIndexX;
    int newIndexY;

    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            particles[i].neighbors.clear();
            // particles[i].neighbors.reserve(14); // not significant

            for (size_t j = 0; j < 9; j++)
            {
                newIndexX = particles[i].k + cellOffset[j][0]; // preallocated cellOffset
                newIndexY = particles[i].l + cellOffset[j][1]; // TODO map grid to cell index in construction to avoid recomputation

                if (newIndexX >= 0 && newIndexY >= 0 && newIndexX < boundingBox[4] && newIndexY < boundingBox[5])
                {
                    uint64_t zIndex = interleaveBits(newIndexX, newIndexY); // less variables

                    for (size_t k = getParticleIndicesZI[zIndex]; k < getParticleIndicesZI[zIndex + 1]; k++) // at() to []
                    {
                        // exit earlier if particles are too far apart, excludes in x-direction half (18) the particles for perfect sampling
                        float dx = particles[i].position.x - particles[k].position.x;
                        if (dx * dx >= h2) continue;

                        float dy = particles[i].position.y - particles[k].position.y;
                        if (dx * dx + dy * dy < h2) {
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
                        float distance = d.x * d.x + d.y * d.y;
                        if (distance < (2.0f * h) * (2.0f * h))
                        {
                            sortedIndicesZI[i].reference->neighbors.push_back(sortedIndicesZI[k].reference);
                        }
                    }
                }
            }
        }
    }
}

void zIndexSortQueryHandleSortImproved(Particle* particles, int numParticles, float h)
{
    // precomputed and allocated
    int h2 = (2.0f * h) * (2.0f * h);
    int newIndexX;
    int newIndexY;
    int cellSizeX = boundingBox[4];
    int cellSizeY = boundingBox[5];
    uint64_t zIndex;
    float dx;
    float dy;
    float mult;

    for (size_t i = 0; i < numParticles; i++)
    {
        // reduce lookups
        float posX = sortedIndicesZI[i].reference->position.x;
        float posY = sortedIndicesZI[i].reference->position.y;
        int cellIdentifierK = sortedIndicesZI[i].reference->k;
        int cellIdentifierL = sortedIndicesZI[i].reference->l;

        if (sortedIndicesZI[i].reference->isFluid)
        {
            sortedIndicesZI[i].reference->neighbors.clear();

            for (size_t j = 0; j < 9; j++)
            {
                newIndexX = cellIdentifierK + cellOffset[j][0]; // preallocted
                newIndexY = cellIdentifierL + cellOffset[j][1];

                if (newIndexX >= 0 && newIndexX < cellSizeX && newIndexY >= 0 && newIndexY < cellSizeY)
                {
                    zIndex = interleaveBits(newIndexX, newIndexY); // less variables // TODO lookuptable for most common?

                    for (size_t k = getParticleIndicesZI[zIndex]; k < getParticleIndicesZI[zIndex + 1]; k++) // at() to []
                    {
                        // exit earlier if particles are too far apart, excludes in x-direction half (18) the particles for perfect sampling
                        dx = posX - sortedIndicesZI[k].reference->position.x;
                        mult = dx * dx;
                        if (mult >= h2) continue;

                        dy = posY - sortedIndicesZI[k].reference->position.y;
                        if (mult + dy * dy < h2) {
                            sortedIndicesZI[i].reference->neighbors.push_back(sortedIndicesZI[k].reference);
                        }
                    }
                }
            }
        }
    }
}

void zIndexSortQueryHandleSortOverCellsImproved(Particle* particles, int numParticles, float h)
{
    // precomputed
    int h2 = (2.0f * h) * (2.0f * h);

    // over cells
    for (size_t i = 0; i < requiredSize - 1; i++)
    {
        int index = getParticleIndicesZI[i];
        if (index >= numParticles) continue;
        int currentK = sortedIndicesZI[index].reference->k;
        int currentL = sortedIndicesZI[index].reference->l;

        // save computing the z cell index a 3/4 of the time
        uint64_t neighborCellsZIndices[] = {
            {interleaveBits(currentK + cellOffset[0][0], currentL + cellOffset[0][1])},
            {interleaveBits(currentK + cellOffset[1][0], currentL + cellOffset[1][1])},
            {interleaveBits(currentK + cellOffset[2][0], currentL + cellOffset[2][1])},
            {interleaveBits(currentK + cellOffset[3][0], currentL + cellOffset[3][1])},
            {interleaveBits(currentK + cellOffset[4][0], currentL + cellOffset[4][1])},
            {interleaveBits(currentK + cellOffset[5][0], currentL + cellOffset[5][1])},
            {interleaveBits(currentK + cellOffset[6][0], currentL + cellOffset[6][1])},
            {interleaveBits(currentK + cellOffset[7][0], currentL + cellOffset[7][1])},
            {interleaveBits(currentK + cellOffset[8][0], currentL + cellOffset[8][1])}
        };

        // over particles in cell
        for (size_t j = index; j < getParticleIndicesZI[i + 1]; j++)
        {
            // Particle* current = sortedIndicesZI[j].reference; // not significant
            if (sortedIndicesZI[j].reference->isFluid)
            {
                sortedIndicesZI[j].reference->neighbors.clear();

                // over neighbor cells
                for (size_t k = 0; k < 9; k++)
                {
                    if (neighborCellsZIndices[k] >= 0 && neighborCellsZIndices[k] < requiredSize) // too many cell for edge cases
                    {
                        // over potential neighbors
                        for (size_t y = getParticleIndicesZI[neighborCellsZIndices[k]]; y < getParticleIndicesZI[neighborCellsZIndices[k] + 1]; y++)
                        {
                            float dx = sortedIndicesZI[j].reference->position.x - sortedIndicesZI[y].reference->position.x;
                            if (dx * dx >= h2) continue;

                            float dy = sortedIndicesZI[j].reference->position.y - sortedIndicesZI[y].reference->position.y;
                            if (dx * dx + dy * dy < h2) {
                                sortedIndicesZI[j].reference->neighbors.push_back(sortedIndicesZI[y].reference);
                            }
                        }
                    }
                }
            }
        }
    }
}