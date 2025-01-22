#include "neighborSearch.h"
#include <iostream>
#include <chrono>

using namespace std::chrono;

vector<vector<Particle*>> cells;
float* boundingBox = new float[6];
vector<int> getParticleIndices;
const int m = 2000;
vector<Particle*> hashTable[m]; // size: 2 * number of particles
int handleArray[m];
vector<vector<Particle*>> compactList;

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

    // add needed cells
    cells.resize(boundingBox[4] * boundingBox[5] + 1);
    for (size_t i = 0; i < boundingBox[4] * boundingBox[5] + 1; i++)
    {
        cells.at(i).clear();
    }
    // compute cell index with (k, l, m)
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
    getParticleIndices.resize(boundingBox[4] * boundingBox[5] + 1);
    fill(getParticleIndices.begin(), getParticleIndices.end(), 0);

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
    // getParticleIndices.at(0) -= 1;
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
                if (cellIndices[j] >= 0 && cellIndices[j] < boundingBox[4] * boundingBox[5]) // valid cell index TODO: change to k,l
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

// Helper function to spread the bits of a 32-bit integer into 64 bits // thx chatgpt
static uint64_t spreadBits(uint32_t x) {
    uint64_t result = x;
    result = (result | (result << 16)) & 0x0000FFFF0000FFFF;
    result = (result | (result << 8)) & 0x00FF00FF00FF00FF;
    result = (result | (result << 4)) & 0x0F0F0F0F0F0F0F0F;
    result = (result | (result << 2)) & 0x3333333333333333;
    result = (result | (result << 1)) & 0x5555555555555555;
    return result;
}

// Function to interleave bits of two 32-bit integers x and y
static uint64_t interleaveBits(uint32_t x, uint32_t y) {
    return (spreadBits(x) | (spreadBits(y) << 1));
}

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
    int requiredSize = biggerXY * biggerXY; // z-index can be alot bigger -> always in 2^n rekursive boxes, unnessesary used memory
    getParticleIndices.resize(requiredSize); // every cell points to the first of its particles
    fill(getParticleIndices.begin(), getParticleIndices.end(), 0);

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
        getParticleIndices.at(particles[i].cellIndex) += 1;
    }
    // accumulate cellIndices
    for (size_t i = 1; i < requiredSize; i++)
    {
        getParticleIndices.at(i) += getParticleIndices.at(i - 1);
    }
    // sort particles with respect to their index with the help of cellIndices
    // insertion sort
    // getParticleIndices.at(0) -= 1; why did i do that?
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
                    for (size_t k = getParticleIndices.at(zIndex); k < getParticleIndices.at(zIndex + 1); k++)
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

int hashFunction(int cellIndexX, int cellIndexY, int sizeHashTable, float d)
{
    return (static_cast<int>(cellIndexX / d * 73856093) xor static_cast<int>(cellIndexY / d * 19349663)) % sizeHashTable;
}

void spatialHashingConstruction(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // reset/ remove old particles
    for (size_t i = 0; i < m; i++)
    {
        hashTable[i].clear();
    }

    // maps grid cell to a hash cell
    for (size_t i = 0; i < numParticles; i++)
    {
        // compute cell index c or cell identifier (k, l, m) for particles
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h)); // conversion fails sometimes on edge cases like 1.000 -> 0
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));

        particles[i].k = k;
        particles[i].l = l;

        // compute hash function i = h(c) or i = h(k, l, m)
        int hashIndex = hashFunction(k, l, m, 1000.f); // / d ? prevent integer overflow
        // store particles in array (hash table) at index i (array of vectors)
        hashTable[hashIndex].push_back(&particles[i]); // store indicies instead of pointer?
    }
}

void spatialHashingQuery(Particle* particles, int numParticles, float h)
{
    // with hash function h(c)
    for (size_t i = 0; i < numParticles; i++) // TODO actually doesnt need this as particles don't change position
    {
        // only fluid needs neighbors
        if (particles[i].isFluid)
        {
            // remove old neighbors
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
                    // cell location in hash table
                    int hashIndex = hashFunction(newIndexX, newIndexY, m, 1000.f);

                    // compare particle distances in the cell with current particle
                    for (size_t k = 0; k < hashTable[hashIndex].size(); k++)
                    {
                        // compute distance
                        Vector2f d = Vector2f(particles[i].position.x - hashTable[hashIndex].at(k)->position.x, particles[i].position.y - hashTable[hashIndex].at(k)->position.y);
                        float distance = sqrt(d.x * d.x + d.y * d.y); // float distanceSquared = d.x * d.x + d.y * d.y; if (distanceSquared < (2.0f * h) * (2.0f * h))
                        
                        // check if neighbor
                        if (distance < 2.0f * h)
                        {
                            // prevent duplicates -> sets?
                            bool duplicate = false;
                            for (size_t z = 0; z < particles[i].neighbors.size(); z++)
                            {
                                if (particles[i].neighbors.at(z)->index == hashTable[hashIndex].at(k)->index)
                                {
                                    duplicate = true;
                                    break;
                                }
                            }
                            if (!duplicate)
                            {
                                // add to neighbors
                                particles[i].neighbors.push_back(hashTable[hashIndex].at(k));
                            }
                        }
                    }
                }
            }
        }
    }
}

void compactHashingConstruction(Particle* particles, int numParticles, float h)
{
    // secondary data structure to store a compact list of non-empty cells
    // hash cells store a handle to corresponding used cell
    // only allocate if used
    // hash table, handle array, compact lists
    // TODO: improvements, preallocate lists with maximally expected number of particles in a cell?, z-sorted secondary structure

    boundingBoxConstruction(particles, numParticles, h);

    // reset/ remove old particles
    fill(handleArray, handleArray + m, 0);
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
        int hashIndex = hashFunction(k, l, m, 1000.f); // / d ? prevent integer overflow

        // store particles in array (hash table) at index i (array of vectors)
        if (handleArray[hashIndex] == 0)
        {
            handleArray[hashIndex] = counter;
            compactList.push_back(vector<Particle*>());
            compactList.at(counter - 1).push_back(&particles[i]);
            counter += 1;
        }
        else
        {
            compactList.at(handleArray[hashIndex] - 1).push_back(&particles[i]);
        }
    }
}

void compactHashingQuery(Particle* particles, int numParticles, float h)
{
    // with hash function h(c)
    for (size_t i = 0; i < numParticles; i++)
    {
        // only fluid needs neighbors
        if (particles[i].isFluid)
        {
            // remove old neighbors
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
                    // cell location in hash table
                    int hashIndex = hashFunction(newIndexX, newIndexY, m, 1000.f);
                    int usedCellIndex = handleArray[hashIndex];

                    // compare particle distances in the cell with current particle
                    for (size_t k = 0; k < compactList.at(usedCellIndex).size(); k++)
                    {
                        // compute distance
                        Vector2f d = Vector2f(particles[i].position.x - compactList.at(usedCellIndex).at(k)->position.x, particles[i].position.y - compactList.at(usedCellIndex).at(k)->position.y);
                        float distance = sqrt(d.x * d.x + d.y * d.y);

                        // check if neighbor
                        if (distance < 2.0f * h)
                        {
                            // prevent duplicates
                            bool duplicate = false;
                            for (size_t z = 0; z < particles[i].neighbors.size(); z++)
                            {
                                if (particles[i].neighbors.at(z)->index == compactList.at(usedCellIndex).at(k)->index)
                                {
                                    duplicate = true;
                                    break;
                                }
                            }
                            if (!duplicate)
                            {
                                // add to neighbors
                                particles[i].neighbors.push_back(compactList.at(usedCellIndex).at(k));
                            }
                        }
                    }
                }
            }
        }
    }
}