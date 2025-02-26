#include "cellLinkedList.h"
#include <iostream>

int lengthCompactList = 0;
int* compactCellIndex = nullptr;
int* compactReference = nullptr;

void cellLinkedListConstruction(Particle* particles, int numParticles, float h)
{
	boundingBoxConstruction(particles, numParticles, h);

    float invCellSize = 1.0f / (2.f * h);

    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);
        particles[i].k = k;
        particles[i].l = l;
        particles[i].cellIndex = interleaveBits(k, l);
    }

    for (size_t i = 1; i < numParticles; i++) // or parallel radix sort, parallel stable merge sort (recommended)
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

    // compact list of all non-empty cells from the sorted particle list
    // for each used cell store cell index and reference to first particle
    // number of particles: difference between particle references
    // 1. compute marker from sorted particle list (1 if cellindex different form preceeding)
    // 2. compute prefix sum (accumulation) parallel scan
    // 3. write particle index to compact cell array at scan location if marker is 1
    // are supposed to be well suited for parallelization

    // parallelization with reduce?
    lengthCompactList = 1;
    for (size_t i = 1; i < numParticles; i++)
    {
        if (particles[i].cellIndex != particles[i - 1].cellIndex) lengthCompactList++;
    }
    if (compactCellIndex != nullptr)
    {
        delete[] compactCellIndex;
        delete[] compactReference;
        compactCellIndex = nullptr;
        compactReference = nullptr;
    }
    compactCellIndex = new int[lengthCompactList];
    compactReference = new int[lengthCompactList + 1];

    int prefixSum = 0; // probably has to be stored for parallelization
    compactCellIndex[0] = particles[0].cellIndex;
    compactReference[0] = 0;

    for (size_t i = 1; i < numParticles; i++)
    {
        if (particles[i].cellIndex != particles[i - 1].cellIndex)
        {
            prefixSum++;
            compactCellIndex[prefixSum] = particles[i].cellIndex;
            compactReference[prefixSum] = i;
        }
    }
    compactReference[prefixSum + 1] = numParticles;
}

void cellLinkedListQueryCountNeighbors(Particle* particles, int numParticles, float h)
{
    if (numNeighbors == nullptr)
    {
        numNeighbors = new int[numParticles];
    }
    memset(&numNeighbors[0], 0, numParticles * sizeof(numNeighbors[0]));

    float h2 = (2.0f * h) * (2.0f * h);
    for (size_t i = 0; i < lengthCompactList; i++)
    {
        for (size_t j = compactReference[i]; j < compactReference[i + 1]; j++)
        {
            if (particles[j].isFluid)
            {
                int neighbors = 0;

                for (size_t h = 0; h < 9; h++)
                {
                    int identifierX = particles[j].k + cellOffset[h][0];
                    int identifierY = particles[j].l + cellOffset[h][1];
                    if (identifierX < boundingBox[4] && identifierY < boundingBox[5])
                    {
                        int cellIndex = interleaveBits(identifierX, identifierY);

                        auto it = lower_bound(compactCellIndex, compactCellIndex + lengthCompactList, cellIndex);
                        if (it == compactCellIndex + lengthCompactList || *it != cellIndex) continue;
                        int reference = it - compactCellIndex;

                        for (size_t k = compactReference[reference]; k < compactReference[reference + 1]; k++)
                        {
                            float dx = particles[j].position.x - particles[k].position.x;
                            float dy = particles[j].position.y - particles[k].position.y;

                            if (dx * dx + dy * dy < h2) {
                                neighbors++;
                            }
                        }
                    }
                }

                numNeighbors[j] = neighbors;
            }
        }
    }
}

void cellLinkedListQuery(Particle* particles, int numParticles, float h) // is probably more efficient in 3D and larger amounts of data
{
    // query over cells
    // cells are searched with with BigMin-LitMax -> other search algorithm, ternary search with fall back to linear search?
    // found neighbors are compressed and stored per particle
    // but uncompress for physics!

    cellLinkedListQueryCountNeighbors(particles, numParticles, h);

    float h2 = (2.0f * h) * (2.0f * h);
    for (size_t i = 0; i < lengthCompactList; i++)
    {
        for (size_t j = compactReference[i]; j < compactReference[i + 1]; j++)
        {
            if (particles[j].isFluid)
            {
                vector<int> neighbors;

                for (size_t h = 0; h < 9; h++)
                {
                    int identifierX = particles[j].k + cellOffset[h][0];
                    int identifierY = particles[j].l + cellOffset[h][1];
                    if (identifierX < boundingBox[4] && identifierY < boundingBox[5])
                    {
                        int cellIndex = interleaveBits(identifierX, identifierY);

                        auto it = lower_bound(compactCellIndex, compactCellIndex + lengthCompactList, cellIndex);
                        if (it == compactCellIndex + lengthCompactList || *it != cellIndex) continue;
                        int reference = it - compactCellIndex;

                        for (size_t k = compactReference[reference]; k < compactReference[reference + 1]; k++)
                        {
                            float dx = particles[j].position.x - particles[k].position.x;
                            float dy = particles[j].position.y - particles[k].position.y;

                            if (dx * dx + dy * dy < h2) {
                                neighbors.push_back(k);
                            }
                        }
                    }
                }
                particles[j].compressedNeighbors = compress(neighbors);
            }
        }
    }
}