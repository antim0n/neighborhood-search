#include <iostream>
#include <unordered_set>
#include <omp.h>
#include "compactHashing.h"

int handleArray[2000];
int hashTableSizeCH = 2000;
int* hashTableCH = nullptr;
vector<vector<int>> compactList;
Handle* sortedIndicesCH = nullptr;
int maxValCH = 10;
int globalCounterCH = 10;


void compactHashingConstruction(Particle* particles, int numParticles, float h) // TODO don't reset complete structure but move changed particles
{
    // secondary data structure to store a compact list of non-empty cells
    // hash cells store a handle to corresponding used cell
    // only allocate if used
    // hash table, handle array, compact lists
    // TODO: improvements, preallocate lists with maximally expected number of particles in a cell?, z-sorted secondary structure

    boundingBoxConstruction(particles, numParticles, h);

    // reset/ remove old particles
    fill(handleArray, handleArray + hashTableSizeCH, 0);
    compactList.clear(); // temporal coherance improvement doesn't need this

    // maps grid cell to a hash cell
    int counter = 1; // to differenatiate from "0" entries
    for (size_t i = 0; i < numParticles; i++)
    {
        // compute cell index c or cell identifier (k, l, m) for particles
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h)); // conversion fails sometimes on edge cases like 1.000 -> 0
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));

        particles[i].k = k;
        particles[i].l = l;

        // compute hash function i = h(c) or i = h(k, l, m)
        int hashIndex = hashFunction(k, l, hashTableSizeCH); // / d ? prevent integer overflow

        // store particles in array (hash table) at index i (array of vectors)
        if (handleArray[hashIndex] == 0)
        {
            handleArray[hashIndex] = counter;
            compactList.push_back(vector<int>());
            compactList.at(counter - 1).push_back(i);
            counter += 1;
        }
        else
        {
            compactList.at(handleArray[hashIndex] - 1).push_back(i);
        }
    }
}

void compactHashingConstructionImproved(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    if (hashTableCH == nullptr)
    {
        hashTableSizeCH = 10 * numParticles; // larger hash table compared to spatial, TODO find good table size
        hashTableCH = new int[hashTableSizeCH];
    }
    memset(&hashTableCH[0], 0, hashTableSizeCH * sizeof(int)); // faster than fill(), TODO most of the table is 0 anyways, can this be improved? maybe also store extra info when it was modified
    
    compactList.clear();
    compactList.reserve(boundingBox[6] / 2); // not significant

    // precompute
    float invCellSize = 1.f / (2.f * h);
    int counter = 1;
    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);

        particles[i].k = k;
        particles[i].l = l;

        int hashIndex = hashFunction(k, l, hashTableSizeCH);

        if (hashTableCH[hashIndex] == 0)
        {
            hashTableCH[hashIndex] = counter;
            compactList.push_back(vector<int>());
            compactList[counter - 1].reserve(8); // preallocate memory to avoid unnecessary reallocations, back() could be used, not significant
            compactList[counter - 1].push_back(i);
            counter += 1;
        }
        else
        {
            compactList[hashTableCH[hashIndex] - 1].push_back(i);
        }
    }
}

void compactHashingConstructionHashCollisionFlagImproved(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    if (hashTableCH == nullptr)
    {
        hashTableSizeCH = 10 * numParticles; // larger hash table compared to spatial, TODO find good table size
        hashTableCH = new int[hashTableSizeCH];
    }
    memset(&hashTableCH[0], 0, hashTableSizeCH * sizeof(int)); // faster than fill(), TODO most of the table is 0 anyways, can this be improved? maybe also store extra info when it was modified

    compactList.clear();
    compactList.reserve(boundingBox[6] / 2); // not significant

    // precompute
    float invCellSize = 1.f / (2.f * h);
    int counter = 1;
    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);

        particles[i].k = k;
        particles[i].l = l;

        int hashIndex = hashFunction(k, l, hashTableSizeCH);

        if (hashTableCH[hashIndex] == 0)
        {
            hashTableCH[hashIndex] = counter;
            compactList.push_back(vector<int>());
            compactList[counter - 1].reserve(8); // preallocate memory to avoid unnecessary reallocations, back() could be used, not significant
            compactList[counter - 1].push_back(i);
            counter += 1;
        }
        else
        {
            int actualIndex = hashTableCH[hashIndex] - 1;
            if ((particles[i].k != particles[compactList[actualIndex].back()].k || particles[i].l != particles[compactList[actualIndex].back()].l) && compactList[actualIndex][0] != 0) // check for hash collisions
            {
                compactList[actualIndex].insert(compactList[actualIndex].begin(), -1);
            }
            compactList[actualIndex].push_back(i);
        }
    }
}

void compactHashingConstructionZSortedImproved(Particle* particles, int numParticles, float h) // TODO hilber curve in comparison
{
    boundingBoxConstruction(particles, numParticles, h);

    if (hashTableCH == nullptr)
    {
        hashTableSizeCH = 10 * numParticles; // larger hash table
        hashTableCH = new int[hashTableSizeCH];
    }
    memset(&hashTableCH[0], 0, hashTableSizeCH * sizeof(int));

    compactList.clear();
    compactList.reserve(boundingBox[6] / 2); // half of max cells as a lot are empty

    // precompute
    float invCellSize = 1.f / (2.f * h);
    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);

        particles[i].k = k;
        particles[i].l = l;
        particles[i].cellIndex = interleaveBits(k, l);
    }

    // insertion sort every step
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

    // insert into hash table
    int counter = 1; // to differenatiate from "0" entries
    for (size_t i = 0; i < numParticles; i++)
    {
        int hashIndex = hashFunction(particles[i].k, particles[i].l, hashTableSizeCH);

        if (hashTableCH[hashIndex] == 0)
        {
            hashTableCH[hashIndex] = counter;
            compactList.push_back(vector<int>());
            compactList[counter - 1].reserve(4); // at() to []
            compactList[counter - 1].push_back(i);
            counter += 1;
        }
        else
        {
            compactList[hashTableCH[hashIndex] - 1].push_back(i);
        }
    }
}

void compactHashingConstructionHandleSort(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // reset/ remove old particles
    fill(handleArray, handleArray + hashTableSizeCH, 0);
    compactList.clear(); // X temporal coherance improvement doesn't need this, temporal coherance not usable
    compactList.reserve(numParticles / 4); // should at least need this amount of cells

    // secondary structure
    if (globalCounterCH == 1)
    {
        if (sortedIndicesCH == nullptr)
        {
            sortedIndicesCH = new Handle[numParticles];
        }
        for (size_t i = 0; i < numParticles; i++)
        {
            sortedIndicesCH[i] = { particles[i].cellIndex, static_cast<int>(i) };
        }
        for (size_t i = 1; i < numParticles; i++)
        {
            Handle current = sortedIndicesCH[i];
            int j = i - 1;
            while (j >= 0 && current.cellIndex < sortedIndicesCH[j].cellIndex)
            {
                sortedIndicesCH[j + 1] = sortedIndicesCH[j];
                j -= 1;
            }
            sortedIndicesCH[j + 1] = current;
        }
        // copy particles twice
        Particle* sortedParticles = new Particle[numParticles];
        for (size_t i = 0; i < numParticles; i++)
        {
            sortedParticles[i] = particles[sortedIndicesCH[i].location];
        }
        copy(sortedParticles, sortedParticles + numParticles, particles);
        delete[] sortedParticles;
        for (size_t i = 0; i < numParticles; i++)
        {
            sortedIndicesCH[i] = { particles[i].cellIndex, static_cast<int>(i) };
        }
    }
    globalCounterCH += 1;
    globalCounterCH %= 100;

    if (globalCounterCH == 1) // TODO: secondary structure
    {
        // z-sorted
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
    globalCounterCH += 1;
    globalCounterCH %= 100;

    // maps grid cell to a hash cell
    int counter = 1; // to differenatiate from "0" entries
    for (size_t i = 0; i < numParticles; i++)
    {
        // compute cell index c or cell identifier (k, l, m) for particles
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h)); // conversion fails sometimes on edge cases like 1.000 -> 0
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));

        particles[i].k = k;
        particles[i].l = l;

        if (globalCounterCH == 1)
        {
            // z-index
            uint32_t x = k;
            uint32_t y = l;
            uint64_t zindex = interleaveBits(x, y);
            particles[i].cellIndex = zindex;
        }

        // compute hash function i = h(c) or i = h(k, l, m)
        int hashIndex = hashFunction(k, l, hashTableSizeCH); // / d ? prevent integer overflow

        // store particles in array (hash table) at index i (array of vectors)
        if (handleArray[hashIndex] == 0)
        {
            handleArray[hashIndex] = counter;
            compactList.push_back(vector<int>());
            compactList.at(counter - 1).reserve(4); // preallocate k entries
            compactList.at(counter - 1).push_back(i);
            counter += 1;
        }
        else
        {
            compactList.at(handleArray[hashIndex] - 1).push_back(i);
        }
    }
}

void compactHashingConstructionHandleSortImproved(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    if (hashTableCH == nullptr)
    {
        hashTableSizeCH = 10 * numParticles; // larger hash table
        hashTableCH = new int[hashTableSizeCH];
    }
    memset(&hashTableCH[0], 0, hashTableSizeCH * sizeof(int));

    compactList.clear();
    compactList.reserve(boundingBox[6] / 2); // half of max cells as a lot are empty

    // precompute
    int counter = 1; // to differenatiate from "0" entries
    float invCellSize = 1.f / (2.f * h);
    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);

        particles[i].k = k;
        particles[i].l = l;
        particles[i].cellIndex = interleaveBits(k, l);

        int hashIndex = hashFunction(k, l, hashTableSizeCH);

        if (hashTableCH[hashIndex] == 0)
        {
            hashTableCH[hashIndex] = counter;
            compactList.push_back(vector<int>());
            compactList[counter - 1].reserve(4); // at() to []
            compactList[counter - 1].push_back(i);
            counter += 1;
        }
        else
        {
            compactList[hashTableCH[hashIndex] - 1].push_back(i);
        }
    }


    if (globalCounterCH == maxValCH)
    {
        // insertion sort every nth step // TODO only for handle structure?
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
        // rebuild handle strucure
        int counter = 0;
        for (size_t i = 0; i < numParticles; i++)
        {
            int hashIndex = hashFunction(particles[i].k, particles[i].l, hashTableSizeCH);

            if (hashTableCH[hashIndex] == 0)
            {
                hashTableCH[hashIndex] = counter;
                compactList.push_back(vector<int>());
                compactList[counter - 1].reserve(4); // at() to []
                compactList[counter - 1].push_back(i);
                counter += 1;
            }
            else
            {
                compactList[hashTableCH[hashIndex] - 1].push_back(i);
            }
        }

        globalCounterCH = 0;
    }

    // build secondary structure // TODO can this be avoided?

    globalCounterCH++;
}

void compactHashingConstructionHandleSortImprovedParallel(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    if (hashTableCH == nullptr)
    {
        hashTableSizeCH = 10 * numParticles;
        hashTableCH = new int[hashTableSizeCH];
    }
    memset(&hashTableCH[0], 0, hashTableSizeCH * sizeof(int));

    compactList.clear();
    compactList.reserve(boundingBox[6] / 2);

    int counter = 1;
    float invCellSize = 1.f / (2.f * h);
    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);

        particles[i].k = k;
        particles[i].l = l;
        particles[i].cellIndex = interleaveBits(k, l);

        int hashIndex = hashFunction(k, l, hashTableSizeCH);

        if (hashTableCH[hashIndex] == 0)
        {
            hashTableCH[hashIndex] = counter;
            compactList.push_back(vector<int>());
            compactList[counter - 1].reserve(4);
            compactList[counter - 1].push_back(i);
            counter += 1;
        }
        else
        {
            compactList[hashTableCH[hashIndex] - 1].push_back(i);
        }
    }

    if (globalCounterCH == maxValCH)
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
        int counter = 0;
        for (size_t i = 0; i < numParticles; i++)
        {
            int hashIndex = hashFunction(particles[i].k, particles[i].l, hashTableSizeCH);

            if (hashTableCH[hashIndex] == 0)
            {
                hashTableCH[hashIndex] = counter;
                compactList.push_back(vector<int>());
                compactList[counter - 1].reserve(4);
                compactList[counter - 1].push_back(i);
                counter += 1;
            }
            else
            {
                compactList[hashTableCH[hashIndex] - 1].push_back(i);
            }
        }
        globalCounterCH = 0;
    }
    globalCounterCH++;
}

void compactHashingQuery(Particle* particles, int numParticles, float h) // TODO Remove not needed function parameters
{
    for (size_t i = 0; i < compactList.size(); i++) // list of used cells is queried
    {
        for (size_t j = 0; j < compactList.at(i).size(); j++)
        {
            int particleIndex = compactList.at(i).at(j);
            if (particles[particleIndex].isFluid)
            {
                particles[particleIndex].neighbors.clear();

                // calculate used cells
                int cellIndices[][2] = { // maybe leave out 0,0 because it is equal to i in the end
                    {-1, 1}, {0, 1}, {1, 1},
                    {-1, 0}, {0, 0}, {1, 0},
                    {-1, -1}, {0, -1}, {1, -1}
                };

                for (size_t k = 0; k < 9; k++)
                {
                    int hashIndex = hashFunction(particles[particleIndex].k + cellIndices[k][0], particles[particleIndex].l + cellIndices[k][1], hashTableSizeCH);
                    int usedCellIndex = handleArray[hashIndex] - 1; // handleArray store Indicies from 1-2000 not 0-1999

                    if (usedCellIndex != -1) // 0 has no entry in the compact List
                    {
                        for (size_t y = 0; y < compactList.at(usedCellIndex).size(); y++)
                        {
                            Vector2f d = particles[particleIndex].position - particles[compactList.at(usedCellIndex).at(y)].position;
                            float distance = d.x * d.x + d.y * d.y;

                            if (distance < (2.0f * h) * (2.0f * h))
                            {
                                bool duplicate = false;
                                for (size_t z = 0; z < particles[particleIndex].neighbors.size(); z++)
                                {
                                    if (particles[particles[particleIndex].neighbors.at(z)].index == particles[compactList.at(usedCellIndex).at(y)].index)
                                    {
                                        duplicate = true;
                                        break;
                                    }
                                }
                                if (!duplicate)
                                {
                                    particles[particleIndex].neighbors.push_back(compactList.at(usedCellIndex).at(y));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void compactHashingQueryImproved(Particle* particles, int numParticles, float h)
{
    float h2 = (2.0f * h) * (2.0f * h); // Precompute

    for (size_t i = 0; i < compactList.size(); i++) // auto is a litte slower
    {
        for (size_t j = 0; j < compactList[i].size(); j++)
        {
            int particleIndex = compactList[i][j];

            if (particles[particleIndex].isFluid)
            {
                particles[particleIndex].neighbors.clear();

                // cell indices have to be computed for every particle separately
                int currentK = particles[particleIndex].k;
                int currentL = particles[particleIndex].l;
                int cellIndices[] =
                {
                    { hashFunction(currentK + cellOffset[0][0], currentL + cellOffset[0][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[1][0], currentL + cellOffset[1][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[2][0], currentL + cellOffset[2][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[3][0], currentL + cellOffset[3][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[4][0], currentL + cellOffset[4][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[5][0], currentL + cellOffset[5][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[6][0], currentL + cellOffset[6][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[7][0], currentL + cellOffset[7][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[8][0], currentL + cellOffset[8][1], hashTableSizeCH) }
                };

                for (size_t k = 0; k < 9; k++)
                {
                    int usedCellIndex = hashTableCH[cellIndices[k]] - 1;

                    if (usedCellIndex != -1)
                    {
                        for (size_t y = 0; y < compactList[usedCellIndex].size(); y++)
                        {
                            int neighborIndex = compactList[usedCellIndex][y];
                            Vector2f d = particles[particleIndex].position - particles[neighborIndex].position;
                            // float dx = compactList[i][j]->position.x - compactList[usedCellIndex][y]->position.x; // slower than vector??
                            // float dy = compactList[i][j]->position.y - compactList[usedCellIndex][y]->position.y;
                            // if (d.x * d.x >= h2) continue; // early exit

                            if (d.x * d.x + d.y * d.y < h2)
                            {
                                bool duplicate = false;
                                for (size_t z = 0; z < particles[particleIndex].neighbors.size(); z++)
                                {
                                    if (particles[particles[particleIndex].neighbors[z]].index == particles[neighborIndex].index)
                                    {
                                        duplicate = true;
                                        break;
                                    }
                                }
                                if (!duplicate)
                                {
                                    particles[particleIndex].neighbors.push_back(neighborIndex);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void compactHashingQueryHashCollisionFlagImproved(Particle * particles, int numParticles, float h)
{
    float h2 = (2.0f * h) * (2.0f * h); // Precompute
    for (size_t i = 0; i < compactList.size(); i++) // auto is a litte slower
    {
        if (compactList[i][0] != -1) // check for hash collision flag, not signifiant improvement
        {
            // precompute neighbor cells, only computed once per cell most of the time
            int currentK = particles[compactList[i][0]].k;
            int currentL = particles[compactList[i][0]].l;
            int cellIndices[] =
            {
                { hashFunction(currentK + cellOffset[0][0], currentL + cellOffset[0][1], hashTableSizeCH) }, // precomputed offset
                { hashFunction(currentK + cellOffset[1][0], currentL + cellOffset[1][1], hashTableSizeCH) },
                { hashFunction(currentK + cellOffset[2][0], currentL + cellOffset[2][1], hashTableSizeCH) },
                { hashFunction(currentK + cellOffset[3][0], currentL + cellOffset[3][1], hashTableSizeCH) },
                { hashFunction(currentK + cellOffset[4][0], currentL + cellOffset[4][1], hashTableSizeCH) },
                { hashFunction(currentK + cellOffset[5][0], currentL + cellOffset[5][1], hashTableSizeCH) },
                { hashFunction(currentK + cellOffset[6][0], currentL + cellOffset[6][1], hashTableSizeCH) },
                { hashFunction(currentK + cellOffset[7][0], currentL + cellOffset[7][1], hashTableSizeCH) },
                { hashFunction(currentK + cellOffset[8][0], currentL + cellOffset[8][1], hashTableSizeCH) }
            };

            for (size_t j = 0; j < compactList[i].size(); j++) // at() seems to be faster than [] here
            {
                // store compactList.at(i).at(j) temporarily, reduces lookup with significant performance boost
                int particleIndex = compactList[i][j];

                if (particles[particleIndex].isFluid) // at to []
                {
                    particles[particleIndex].neighbors.clear();

                    for (size_t k = 0; k < 9; k++)
                    {
                        int usedCellIndex = hashTableCH[cellIndices[k]] - 1;

                        if (usedCellIndex != -1)
                        {
                            for (size_t y = 0; y < compactList[usedCellIndex].size(); y++)
                            {
                                int neighborIndex = compactList[usedCellIndex][y];
                                if (neighborIndex == -1) continue; // skip hash collision flags

                                Vector2f d = particles[particleIndex].position - particles[neighborIndex].position;

                                if (d.x * d.x + d.y * d.y < h2)
                                {
                                    bool duplicate = false;
                                    for (size_t z = 0; z < particles[particleIndex].neighbors.size(); z++)
                                    {
                                        if (particles[particles[particleIndex].neighbors[z]].index == particles[neighborIndex].index)
                                        {
                                            duplicate = true;
                                            break;
                                        }
                                    }
                                    if (!duplicate)
                                    {
                                        particles[particleIndex].neighbors.push_back(neighborIndex);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else // hash collision in cell occured
        {
            for (size_t j = 1; j < compactList[i].size(); j++) // skip nullpointer
            {
                int particleIndex = compactList[i][j];
                if (particleIndex == -1) continue; // skip hash collision flags

                if (particles[particleIndex].isFluid)
                {
                    particles[particleIndex].neighbors.clear();

                    // cell indices have to be computed for every particle separately
                    int currentK = particles[particleIndex].k;
                    int currentL = particles[particleIndex].l;
                    int cellIndices[] =
                    {
                        { hashFunction(currentK + cellOffset[0][0], currentL + cellOffset[0][1], hashTableSizeCH) },
                        { hashFunction(currentK + cellOffset[1][0], currentL + cellOffset[1][1], hashTableSizeCH) },
                        { hashFunction(currentK + cellOffset[2][0], currentL + cellOffset[2][1], hashTableSizeCH) },
                        { hashFunction(currentK + cellOffset[3][0], currentL + cellOffset[3][1], hashTableSizeCH) },
                        { hashFunction(currentK + cellOffset[4][0], currentL + cellOffset[4][1], hashTableSizeCH) },
                        { hashFunction(currentK + cellOffset[5][0], currentL + cellOffset[5][1], hashTableSizeCH) },
                        { hashFunction(currentK + cellOffset[6][0], currentL + cellOffset[6][1], hashTableSizeCH) },
                        { hashFunction(currentK + cellOffset[7][0], currentL + cellOffset[7][1], hashTableSizeCH) },
                        { hashFunction(currentK + cellOffset[8][0], currentL + cellOffset[8][1], hashTableSizeCH) }
                    };

                    for (size_t k = 0; k < 9; k++)
                    {
                        int usedCellIndex = hashTableCH[cellIndices[k]] - 1;

                        if (usedCellIndex != -1)
                        {
                            for (size_t y = 0; y < compactList[usedCellIndex].size(); y++)
                            {
                                int neighborIndex = compactList[usedCellIndex][y];
                                if (neighborIndex == -1) continue; // skip hash collision flags

                                Vector2f d = particles[particleIndex].position - particles[neighborIndex].position;
                                if (d.x * d.x + d.y * d.y < h2)
                                {
                                    bool duplicate = false;
                                    for (size_t z = 0; z < particles[particleIndex].neighbors.size(); z++)
                                    {
                                        if (particles[particles[particleIndex].neighbors[z]].index == particles[neighborIndex].index)
                                        {
                                            duplicate = true;
                                            break;
                                        }
                                    }
                                    if (!duplicate)
                                    {
                                        particles[particleIndex].neighbors.push_back(neighborIndex);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void compactHashingQueryCountNeighbors(Particle* particles, int numParticles, float h)
{
    if (numNeighbors == nullptr)
    {
        numNeighbors = new int[numParticles];
    }
    memset(&numNeighbors[0], 0, numParticles * sizeof(int));

    float h2 = (2.0f * h) * (2.0f * h);

    // #pragma omp parallel for
    for (int i = 0; i < compactList.size(); i++)
    {
        for (size_t j = 0; j < compactList[i].size(); j++)
        {
            int particleIndex = compactList[i][j];

            if (particles[particleIndex].isFluid)
            {
                int currentK = particles[particleIndex].k;
                int currentL = particles[particleIndex].l;
                int cellIndices[] =
                {
                    { hashFunction(currentK + cellOffset[0][0], currentL + cellOffset[0][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[1][0], currentL + cellOffset[1][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[2][0], currentL + cellOffset[2][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[3][0], currentL + cellOffset[3][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[4][0], currentL + cellOffset[4][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[5][0], currentL + cellOffset[5][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[6][0], currentL + cellOffset[6][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[7][0], currentL + cellOffset[7][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[8][0], currentL + cellOffset[8][1], hashTableSizeCH) }
                };

                unordered_set<int> neighbors;
                neighbors.reserve(20);

                for (size_t k = 0; k < 9; k++)
                {
                    int usedCellIndex = hashTableCH[cellIndices[k]] - 1;

                    if (usedCellIndex != -1)
                    {
                        for (size_t y = 0; y < compactList[usedCellIndex].size(); y++)
                        {
                            int neighborIndex = compactList[usedCellIndex][y];
                            float dx = particles[particleIndex].position.x - particles[neighborIndex].position.x;
                            float dy = particles[particleIndex].position.y - particles[neighborIndex].position.y;

                            if (dx * dx + dy * dy < h2)
                            {
                                neighbors.insert(neighborIndex);
                            }
                        }
                    }
                }
                numNeighbors[particleIndex] = neighbors.size();
            }
        }
    }
}

void compactHashingQueryImprovedParallel(Particle* particles, int numParticles, float h)
{
    compactHashingQueryCountNeighbors(particles, numParticles, h);

    float h2 = (2.0f * h) * (2.0f * h);

    // #pragma omp parallel for
    for (int i = 0; i < compactList.size(); i++)
    {
        for (size_t j = 0; j < compactList[i].size(); j++)
        {
            int particleIndex = compactList[i][j];

            if (particles[particleIndex].isFluid)
            {
                particles[particleIndex].neighbors.clear();
                particles[particleIndex].neighbors.resize(numNeighbors[particleIndex]);
                int counter = 0;

                int currentK = particles[particleIndex].k;
                int currentL = particles[particleIndex].l;
                int cellIndices[] =
                {
                    { hashFunction(currentK + cellOffset[0][0], currentL + cellOffset[0][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[1][0], currentL + cellOffset[1][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[2][0], currentL + cellOffset[2][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[3][0], currentL + cellOffset[3][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[4][0], currentL + cellOffset[4][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[5][0], currentL + cellOffset[5][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[6][0], currentL + cellOffset[6][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[7][0], currentL + cellOffset[7][1], hashTableSizeCH) },
                    { hashFunction(currentK + cellOffset[8][0], currentL + cellOffset[8][1], hashTableSizeCH) }
                };

                for (size_t k = 0; k < 9; k++)
                {
                    int usedCellIndex = hashTableCH[cellIndices[k]] - 1;

                    if (usedCellIndex != -1)
                    {
                        for (size_t y = 0; y < compactList[usedCellIndex].size(); y++)
                        {
                            int neighborIndex = compactList[usedCellIndex][y];
                            Vector2f d = particles[particleIndex].position - particles[neighborIndex].position;

                            if (d.x * d.x + d.y * d.y < h2)
                            {
                                bool duplicate = false;
                                for (size_t z = 0; z < particles[particleIndex].neighbors.size(); z++)
                                {
                                    if (particles[particles[particleIndex].neighbors[z]].index == particles[neighborIndex].index)
                                    {
                                        duplicate = true;
                                        break;
                                    }
                                }
                                if (!duplicate)
                                {
                                    particles[particleIndex].neighbors[counter] = neighborIndex;
                                    counter++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}