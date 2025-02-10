#include "compactHashing.h"

const int m = 2000; // should change when number of fluid particles are changed
int handleArray[m];
vector<vector<Particle*>> compactList;
Handle* sortedIndicesCH = nullptr;
int maxValCH = 10;
int globalCounterCH = 100;


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
        int hashIndex = hashFunction(k, l, m); // / d ? prevent integer overflow

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
void compactHashingConstructionZSorted(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // reset/ remove old particles
    fill(handleArray, handleArray + m, 0);
    compactList.clear(); // X temporal coherance improvement doesn't need this, temporal coherance not usable
    compactList.reserve(numParticles / 4); // should at least need this amount of cells

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
        int hashIndex = hashFunction(k, l, m); // / d ? prevent integer overflow

        // store particles in array (hash table) at index i (array of vectors)
        if (handleArray[hashIndex] == 0)
        {
            handleArray[hashIndex] = counter;
            compactList.push_back(vector<Particle*>());
            compactList.at(counter - 1).reserve(4); // preallocate k entries
            compactList.at(counter - 1).push_back(&particles[i]);
            counter += 1;
        }
        else
        {
            compactList.at(handleArray[hashIndex] - 1).push_back(&particles[i]);
        }
    }
}

void compactHashingConstructionHandleSort(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // reset/ remove old particles
    fill(handleArray, handleArray + m, 0);
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
            sortedIndicesCH[i] = { particles[i].cellIndex, &particles[i] };
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
            sortedParticles[i] = *sortedIndicesCH[i].reference;
        }
        copy(sortedParticles, sortedParticles + numParticles, particles);
        delete[] sortedParticles;
        for (size_t i = 0; i < numParticles; i++)
        {
            sortedIndicesCH[i] = { particles[i].cellIndex, &particles[i] };
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
        int hashIndex = hashFunction(k, l, m); // / d ? prevent integer overflow

        // store particles in array (hash table) at index i (array of vectors)
        if (handleArray[hashIndex] == 0)
        {
            handleArray[hashIndex] = counter;
            compactList.push_back(vector<Particle*>());
            compactList.at(counter - 1).reserve(4); // preallocate k entries
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
    for (size_t i = 0; i < compactList.size(); i++)
    {
        for (size_t j = 0; j < compactList.at(i).size(); j++)
        {
            if (compactList.at(i).at(j)->isFluid)
            {
                compactList.at(i).at(j)->neighbors.clear();

                // calculate used cells
                int cellIndices[][2] = { // maybe leave out 0,0 because it is equal to i in the end
                    {-1, 1}, {0, 1}, {1, 1},
                    {-1, 0}, {0, 0}, {1, 0},
                    {-1, -1}, {0, -1}, {1, -1}
                };

                for (size_t k = 0; k < 9; k++)
                {
                    int hashIndex = hashFunction(compactList.at(i).at(j)->k + cellIndices[k][0], compactList.at(i).at(j)->l + cellIndices[k][1], m);
                    int usedCellIndex = handleArray[hashIndex] - 1; // handleArray store Indicies from 1-2000 not 0-1999

                    if (usedCellIndex != -1) // 0 has no entry in the compact List
                    {
                        for (size_t y = 0; y < compactList.at(usedCellIndex).size(); y++)
                        {
                            Vector2f d = compactList.at(i).at(j)->position - compactList.at(usedCellIndex).at(y)->position;
                            float distance = sqrt(d.x * d.x + d.y * d.y);

                            if (distance < 2.0f * h)
                            {
                                bool duplicate = false;
                                for (size_t z = 0; z < compactList.at(i).at(j)->neighbors.size(); z++)
                                {
                                    if (compactList.at(i).at(j)->neighbors.at(z)->index == compactList.at(usedCellIndex).at(y)->index)
                                    {
                                        duplicate = true;
                                        break;
                                    }
                                }
                                if (!duplicate)
                                {
                                    compactList.at(i).at(j)->neighbors.push_back(compactList.at(usedCellIndex).at(y));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}