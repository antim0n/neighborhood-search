#include <iostream>
#include <omp.h>
#include "spatialHashing.h"

int hashTableSizeSH = 0;
vector<Particle*>* hashTableSH = nullptr; // size: 2 * number of particles

void spatialHashingConstruction(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);
    // allocate memory once
    if (hashTableSH == nullptr)
    {
        hashTableSizeSH = 2 * numParticles;
        hashTableSH = new vector<Particle*>[hashTableSizeSH];
    }

    // reset/ remove old particles
    for (size_t i = 0; i < hashTableSizeSH; i++)
    {
        hashTableSH[i].clear();
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
        int hashIndex = hashFunction(k, l, hashTableSizeSH); // prevent integer overflow
        // store particles in array (hash table) at index i (array of vectors)
        hashTableSH[hashIndex].push_back(&particles[i]); // store indicies instead of pointer?
    }
}

void spatialHashingConstructionImproved(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);
    // allocate memory once
    if (hashTableSH == nullptr)
    {
        //// https://cseweb.ucsd.edu/~kube/cls/100/Lectures/lec16/lec16-8.html
        //// https://www.geeksforgeeks.org/c-program-to-check-prime-number/
        //// https://dl.acm.org/doi/abs/10.1145/356643.356645
        //int n = static_cast<int>(2.0f * numParticles); // bigger hash tables than 2 * numParticles do not seem to increase performance much
        //int j = n;
        //while(true)
        //{
        //    int cnt = 0;
        //    for (int i = 1; i <= j; i++)
        //    {
        //        if (j % i == 0)
        //            cnt++;
        //    }
        //    if (cnt > 2)
        //    {
        //        j++;
        //        continue;
        //    }
        //    else
        //    {
        //        hashTableSize = j;
        //        break;
        //    }
        //}
        // hashTableSize = pow(2, static_cast<int>(log2(2 * numParticles) + 1)); // power of 2 chat gpt seems to prefer this for a "better modulo performance"
        // hashTableSize = 200000; // power of 10
        hashTableSizeSH = 2 * numParticles;
        hashTableSH = new vector<Particle*>[hashTableSizeSH];
    }

    for (size_t i = 0; i < hashTableSizeSH; i++)
    {
        hashTableSH[i].clear(); // TODO clear does set the capacity to 0!!
        // hashTableSH[i].reserve(4); // not a significant improvement as too much empty cells, too much slows down the query
    }

    // precompute and allocate
    float invCellSize = 1.0f / (2.f * h); // TODO maybe  remove from function, with initializer
    int k;
    int l;

    for (size_t i = 0; i < numParticles; i++)
    {
        k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);

        particles[i].k = k;
        particles[i].l = l;

        int hashIndex = hashFunction(k, l, hashTableSizeSH);
        hashTableSH[hashIndex].push_back(&particles[i]); // store indicies instead of pointer? is a lot slower
    }
}

void spatialHashingQuery(Particle* particles, int numFluidParticles, float h)
{
    // with hash function h(c)
    for (size_t i = 0; i < numFluidParticles; i++)
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
                int hashIndex = hashFunction(newIndexX, newIndexY, hashTableSizeSH);

                // compare particle distances in the cell with current particle
                for (size_t k = 0; k < hashTableSH[hashIndex].size(); k++)
                {
                    // compute distance
                    Vector2f d = Vector2f(particles[i].position.x - hashTableSH[hashIndex].at(k)->position.x, particles[i].position.y - hashTableSH[hashIndex].at(k)->position.y);
                    float distance = d.x * d.x + d.y * d.y;
                    // check if neighbor
                    if (distance < (2.0f * h) * (2.0f * h))
                    {
                        // prevent duplicates -> sets?
                        bool duplicate = false;
                        for (size_t z = 0; z < particles[i].neighbors.size(); z++)
                        {
                            if (particles[i].neighbors.at(z)->index == hashTableSH[hashIndex].at(k)->index)
                            {
                                duplicate = true;
                                break;
                            }
                        }
                        if (!duplicate)
                        {
                            // add to neighbors
                            particles[i].neighbors.push_back(hashTableSH[hashIndex].at(k));
                        }
                    }
                }
            }
        }
    }
}

void spatialHashingQueryImproved(Particle* particles, int numFluidParticles, float h)
{
    // precompute and allocate
    float h2 = (2.0f * h) * (2.0f * h);
    int newIndexX;
    int newIndexY;
    int hashIndex;
    float dx;
    float dy;

    for (size_t i = 0; i < numFluidParticles; i++)
    {
        particles[i].neighbors.clear();
        // particles[i].neighbors.reserve(20); // again, does not do much other than improving the first iteration, and slightly slowing the rest

        for (size_t j = 0; j < 9; j++)
        {
            newIndexX = particles[i].k + cellOffset[j][0];
            newIndexY = particles[i].l + cellOffset[j][1];

            if (newIndexX >= 0 && newIndexX < boundingBox[4] && newIndexY >= 0 && newIndexY < boundingBox[5])
            {
                hashIndex = hashFunction(newIndexX, newIndexY, hashTableSizeSH);

                for (size_t k = 0; k < hashTableSH[hashIndex].size(); k++)
                {
                    dx = particles[i].position.x - hashTableSH[hashIndex][k]->position.x;
                    // if (dx * dx >= h2) continue; // does not do much, removing the vector was what helped
                    dy = particles[i].position.y - hashTableSH[hashIndex][k]->position.y;

                    if (dx * dx + dy * dy < h2)
                    {
                        // using an unordered set to check for duplicates slowed it down a lot
                        /*if (visitedNeighbors.insert(hashTableSH[hashIndex][k]->index).second)
                        {
                            particles[i].neighbors.push_back(hashTableSH[hashIndex][k]);
                        }
                        visitedNeighbors.clear();*/

                        bool duplicate = false;
                        for (size_t z = 0; z < particles[i].neighbors.size(); z++)
                        {
                            if (particles[i].neighbors[z]->index == hashTableSH[hashIndex][k]->index)
                            {
                                duplicate = true;
                                break;
                            }
                        }
                        if (!duplicate)
                        {
                            particles[i].neighbors.push_back(hashTableSH[hashIndex][k]);
                        }
                    }
                }
            }
        }
    }
}

void spatialHashingQueryOverCellsImproved(Particle* particles, int numFluidParticles, float h)
{
    // seems to be "impossible" as of now
    // due to hash collisions and for now not knowing if one has occured in the current cell, a single set of 9 neighbor cells can not be computed
    // an arbitrary number of particles from many cells could be in one hash cell
    // the cellindex would have to be stored and as soon as it differs the neighbor cells have to be recomputed
    // also as the hash table has many empty cells, iterating over them takes unnecessary time
}