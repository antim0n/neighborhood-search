#include "spatialHashing.h"

const int m = 2000; // should change when number of fluid particles are changed
vector<Particle*> hashTable[m]; // size: 2 * number of particles

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
        int hashIndex = hashFunction(k, l, m); // prevent integer overflow
        // store particles in array (hash table) at index i (array of vectors)
        hashTable[hashIndex].push_back(&particles[i]); // store indicies instead of pointer?
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
                int hashIndex = hashFunction(newIndexX, newIndexY, m);

                // compare particle distances in the cell with current particle
                for (size_t k = 0; k < hashTable[hashIndex].size(); k++)
                {
                    // compute distance
                    Vector2f d = Vector2f(particles[i].position.x - hashTable[hashIndex].at(k)->position.x, particles[i].position.y - hashTable[hashIndex].at(k)->position.y);
                    float distance = d.x * d.x + d.y * d.y;
                    // check if neighbor
                    if (distance < (2.0f * h) * (2.0f * h))
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