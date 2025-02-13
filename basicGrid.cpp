#include "basicGrid.h"

vector<vector<Particle*>> cells;

void gridConstruction(Particle* particles, int numParticles, float h) // looks pretty optimal
{
    // uniform grid
    // particles are stored in a cell, without sorting
    boundingBoxConstruction(particles, numParticles, h);

    // add needed cells
    cells.resize(boundingBox[4] * boundingBox[5] + 1);
    for (size_t i = 0; i < boundingBox[4] * boundingBox[5] + 1; i++)
    {
        cells.at(i).clear();
        // cells.at(i).reserve(8); // reserving should not improve performance as many cells are empty does
    }

    // compute cell index with (k, l, m)
    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h));
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));
        // particles[i].k = k;
        // particles[i].l = l;
        particles[i].cellIndex = k + l * boundingBox[4];
        cells.at(particles[i].cellIndex).push_back(&particles[i]); // increased runtime in first iteration due to one-time allocations
    }
}

void gridQuery(Particle* particles, int numFluidParticles, float h)
{
    for (size_t i = 0; i < numFluidParticles; i++)
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

void gridQueryGenerallyImproved(Particle* particles, int numFluidParticles, float h)
{
    for (size_t i = 0; i < numFluidParticles; i++)
    {
        particles[i].neighbors.clear();
        // reserve only makes sense for the first iteration, once a particle had neighbors the memory is sufficient
        // but the function does nothing if there is already enough space
        // too much memory might decrease cache hit rates -> 13 is max number of neighbors in perfect sampling
        particles[i].neighbors.reserve(15);

        // compute cell indices
        // could not think of a reason using k,l for the cell index calculation woulds be faster TODO
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
            // only caps beginning and end, could be left out of boundary is always one cell bigger than fluid TODO
            if (cellIndices[j] >= 0 && cellIndices[j] < boundingBox[4] * boundingBox[5])
            {
                for (size_t k = 0; k < cells.at(cellIndices[j]).size(); k++)
                {
                    Vector2f d = particles[i].position - cells.at(cellIndices[j]).at(k)->position;
                    // not using sqrt() brings pretty significant improvement of approx. 20ms for 100k particles
                    float distance = d.x * d.x + d.y * d.y;
                    if (distance < (2.0f * h) * (2.0f * h))
                    {
                        particles[i].neighbors.push_back(cells.at(cellIndices[j]).at(k)); // increased runtime in first iteration due to one-time allocations
                    }
                }
            }
        }
    }
}

void gridQueryOverCells(Particle* particles, int numFluidParticles, float h)
{
    // a lot of empty cells but less cell index calculation TODO
}