#include "basicGrid.h"

vector<vector<Particle*>> cells;

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
        // particles[i].k = k;
        // particles[i].l = l;
        particles[i].cellIndex = k + l * boundingBox[4];
        cells.at(particles[i].cellIndex).push_back(&particles[i]); // increased runtime in first iteration due to one-time allocations
    }
}

void gridConstructionImproved(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // precompute
    int numberOfCells = boundingBox[6];
    // only resize when nessecary
    if (cells.size() != numberOfCells) {
        cells.resize(numberOfCells);
    }
    for (size_t i = 0; i < numberOfCells; i++)
    {
        cells[i].clear();
        // cells.at(i).reserve(4); // reserving should not improve performance as many cells are empty does
    }

    // precompute
    float invCellSize = 1.0f / (2.f * h);
    float offsetX = boundingBox[0];
    float offsetY = boundingBox[2];

    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - offsetX) * invCellSize);
        int l = static_cast<int>((particles[i].position.y - offsetY) * invCellSize);
        particles[i].cellIndex = k + l * boundingBox[4];
        cells[particles[i].cellIndex].push_back(&particles[i]);
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

void gridQueryOverCells(float h)
{
    // a lot of empty cells but less cell index calculation
    // spatial locality probably worse
    for (size_t i = 0; i < boundingBox[4] * boundingBox[5] + 1; i++)
    {
        if (cells.at(i).empty())
        {
            continue;
        }

        int index = cells.at(i).at(0)->cellIndex;
        int cellIndices[] = { index,
            index + 1,
            index - 1,
            index + boundingBox[4],
            index - boundingBox[4],
            index + boundingBox[4] + 1,
            index + boundingBox[4] - 1,
            index - boundingBox[4] + 1,
            index - boundingBox[4] - 1
        };

        for (size_t j = 0; j < cells.at(i).size(); j++)
        {
            cells.at(i).at(j)->neighbors.clear();

            for (size_t k = 0; k < 9; k++)
            {
                if (cellIndices[k] >= 0 && cellIndices[k] < boundingBox[4] * boundingBox[5])
                {
                    for (size_t y = 0; y < cells.at(cellIndices[k]).size(); y++)
                    {
                        Vector2f d = cells.at(i).at(j)->position - cells.at(cellIndices[k]).at(y)->position;
                        float distance = d.x * d.x + d.y * d.y;
                        if (distance < (2.0f * h) * (2.0f * h))
                        {
                            cells.at(i).at(j)->neighbors.push_back(cells.at(cellIndices[k]).at(y));
                        }
                    }
                }
            }
        }
    }
}

void gridQueryImproved(Particle* particles, int numFluidParticles, float h)
{
    // precompute (not significant)
    int bboxX = boundingBox[4];
    int bboxY = boundingBox[5];
    float h2 = (2.0f * h) * (2.0f * h);

    for (size_t i = 0; i < numFluidParticles; i++)
    {
        particles[i].neighbors.clear();
        // reserve only makes sense for the first iteration, once a particle had neighbors the memory is sufficient
        // but the function does nothing if there is already enough space
        // too much memory might decrease cache hit rates -> 13 is max number of neighbors in perfect sampling
        particles[i].neighbors.reserve(14);

        // compute cell indices
        // could not think of a reason using k,l for the cell index calculation woulds be faster
        
        // precompute (not significant)
        int index = particles[i].cellIndex;
        int cellIndices[] = { particles[i].cellIndex,
            index + 1,
            index - 1,
            index + bboxX,
            index - bboxX,
            index + bboxX + 1,
            index + bboxX - 1,
            index - bboxX + 1,
            index - bboxX - 1
        };

        // check particles in all adjacent cells
        for (size_t j = 0; j < 9; j++)
        {
            // only caps beginning and end, could be left out of boundary is always one cell bigger than fluid
            // does lead to checking unnessecary cells if fluid is on the grid border
            // at the moment not relevant as scene has boundary all around
            if (cellIndices[j] >= 0 && cellIndices[j] < bboxX * bboxY)
            {
                // changed acces with .at() to direct indexing
                for (size_t k = 0; k < cells[cellIndices[j]].size(); k++)
                {
                    float dx = particles[i].position.x - cells[cellIndices[j]][k]->position.x;
                    if (dx * dx >= h2) continue;  // Skip unnecessary checks if x-component is too large (not significant)

                    float dy = particles[i].position.y - cells[cellIndices[j]][k]->position.y;
                    float distance = dx * dx + dy * dy;
                    if (distance < h2) {
                        particles[i].neighbors.push_back(cells[cellIndices[j]][k]);
                    }
                }
            }
        }
    }
}