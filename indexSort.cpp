#include <iostream>
#include <omp.h>
#include "indexSort.h"

vector<int> getParticleIndicesI;
Particle* sortedParticlesI = nullptr;

void indexSortConstructionCountingSort(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    if (sortedParticlesI == nullptr)
    {
        sortedParticlesI = new Particle[numParticles];
    }

    getParticleIndicesI.resize(boundingBox[4] * boundingBox[5] + 1);
    fill(getParticleIndicesI.begin(), getParticleIndicesI.end(), 0);

    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h));
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));
        particles[i].cellIndex = k + l * boundingBox[4];
        getParticleIndicesI.at(particles[i].cellIndex) += 1;
    }
    for (size_t i = 1; i < boundingBox[4] * boundingBox[5] + 1; i++)
    {
        getParticleIndicesI.at(i) += getParticleIndicesI.at(i - 1);
    }

    // copy particles to the correct position in a new array and copy sorted array back to particles
    // counting sort (slightly faster than insertion sort)
    for (size_t i = 0; i < numParticles; i++)
    {
        // -- for the correct query index structure
        sortedParticlesI[--getParticleIndicesI.at(particles[i].cellIndex)] = particles[i];
    }
    copy(sortedParticlesI, sortedParticlesI + numParticles, particles);
}

int getMax(Particle* arr, int size)
{
    int m = arr[0].cellIndex;
    for (int i = 1; i < size; i++)
    {
        getParticleIndicesI.at(arr[i].cellIndex) -= 1;
        if (arr[i].cellIndex > m)
        {
            m = arr[i].cellIndex;
        }
    }
    return m;
}

void countingSort(Particle* arr, int size, int exp)
{
    vector<Particle> output(size);
    int count[10] = { 0 };

    for (int i = 0; i < size; i++) {
        int digit = static_cast<int>(arr[i].cellIndex / exp) % 10;
        count[digit]++;
    }

    for (int i = 1; i < 10; i++) {
        count[i] += count[i - 1];
    }

    for (int i = size - 1; i >= 0; i--) {
        int digit = static_cast<int>(arr[i].cellIndex / exp) % 10;
        output[count[digit] - 1] = arr[i];
        count[digit]--;
    }

    for (int i = 0; i < size; i++) {
        arr[i] = output[i];
    }
}

void radixSort(Particle* arr, int size)
{
    int maxCellIndex = getMax(arr, size);

    for (int exp = 1; maxCellIndex / exp > 0; exp *= 10)
    {
        countingSort(arr, size, exp);
    }
}

void merge(Particle* arr, int left, int mid, int right)
{
    int n1 = mid - left + 1;
    int n2 = right - mid;

    vector<Particle> L(n1), R(n2);

    for (int i = 0; i < n1; i++)
        L[i] = arr[left + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[mid + 1 + j];

    int i = 0, j = 0;
    int k = left;

    while (i < n1 && j < n2) {
        if (L[i].cellIndex <= R[j].cellIndex) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}
void mergeSort(Particle* arr, int left, int right)
{
    if (left >= right)
    {
        return;
    }

    int mid = left + (right - left) / 2;
    mergeSort(arr, left, mid);
    mergeSort(arr, mid + 1, right);
    merge(arr, left, mid, right);
}

void insertionSort(Particle* arr, int size)
{
    for (size_t i = 1; i < size; i++)
    {
        Particle current = arr[i];
        getParticleIndicesI.at(current.cellIndex) -= 1;
        int j = i - 1;
        while (j >= 0 && current.cellIndex < arr[j].cellIndex)
        {
            arr[j + 1] = arr[j];
            j -= 1;
        }
        arr[j + 1] = current;
    }
}

// index sort takes particles as input, sorts them and returns an integer array as output
void indexSortConstructionCompareSorting(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    // compute cell index with (k, l, m)
    getParticleIndicesI.resize(boundingBox[4] * boundingBox[5] + 1);
    fill(getParticleIndicesI.begin(), getParticleIndicesI.end(), 0);

    for (size_t i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) / (2.f * h)); // conversion fails sometimes on edge cases like 1.000 -> 0
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) / (2.f * h));
        particles[i].cellIndex = k + l * boundingBox[4];
        // increment cellIndices
        getParticleIndicesI.at(particles[i].cellIndex) += 1;
    }
    // accumulate cellIndices
    for (size_t i = 1; i < boundingBox[4] * boundingBox[5] + 1; i++)
    {
        getParticleIndicesI.at(i) += getParticleIndicesI.at(i - 1);
    }

    // sort particles with respect to their index with the help of cellIndices
    // insertion sort pretty fast if almost sorted (very slow first iteration)
    insertionSort(particles, numParticles);

    // radix sort with counting sort https://www.geeksforgeeks.org/radix-sort/ (is significantly slower as insertion sort if particles are almost sorted)
    // radixSort(particles, numParticles);

    // merge sort https://www.geeksforgeeks.org/merge-sort/ (even slower than radix sort)
    /*for (size_t i = 0; i < numParticles; i++)
    {
        getParticleIndicesI.at(particles[i].cellIndex) -= 1;
    }
    mergeSort(particles, 0, numParticles - 1);*/
}

void indexSortConstructionInsertionSortImproved(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    getParticleIndicesI.resize(boundingBox[6] + 1); // boundingBox[6] contains number of cells
    memset(&getParticleIndicesI[0], 0, (boundingBox[6] + 1) * sizeof(getParticleIndicesI[0])); // memset faster than fill()

    // getParticleIndicesI = vector<int>(boundingBox[6] + 1, 0); // instead of resize and fill (slightly worse)

    // precompute
    float invCellSize = 1.0f / (2.f * h);
    int k;
    int l;
    for (size_t i = 0; i < numParticles; i++)
    {
        k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);
        particles[i].cellIndex = k + l * boundingBox[4];
        getParticleIndicesI[particles[i].cellIndex] += 1; // [] instead of at()
    }

    for (size_t i = 1; i < boundingBox[6] + 1; i++)
    {
        getParticleIndicesI[i] += getParticleIndicesI[i - 1];
    }

    // slower than counting sort :( might get better with the secondary structure where only handles are sorted (less copying for insertion sort)
    for (size_t i = 1; i < numParticles; i++)
    {
        Particle current = particles[i];
        getParticleIndicesI[current.cellIndex] -= 1;
        int j = i - 1;
        while (j >= 0 && current.cellIndex < particles[j].cellIndex)
        {
            particles[j + 1] = particles[j];
            j -= 1;
        }
        particles[j + 1] = current;
    }
}

void indexSortConstructionCountingSortImproved(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    getParticleIndicesI.resize(boundingBox[6] + 1); // boundingBox[6] contains number of cells
    memset(&getParticleIndicesI[0], 0, (boundingBox[6] + 1) * sizeof(getParticleIndicesI[0])); // memset faster than fill()

    // getParticleIndicesI = vector<int>(boundingBox[6] + 1, 0); // instead of resize and fill (slightly worse)

    // precompute
    float invCellSize = 1.0f / (2.f * h);
    int k;
    int l;
    for (size_t i = 0; i < numParticles; i++)
    {
        k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);
        particles[i].k = k;
        particles[i].l = l;
        particles[i].cellIndex = k + l * boundingBox[4];
        getParticleIndicesI[particles[i].cellIndex] += 1; // [] instead of at()
    }

    for (size_t i = 1; i < boundingBox[6] + 1; i++)
    {
        getParticleIndicesI[i] += getParticleIndicesI[i - 1];
    }

    // countingSort
    if (sortedParticlesI == nullptr)
    {
        sortedParticlesI = new Particle[numParticles];
    }
    for (size_t i = 0; i < numParticles; i++)
    {
        sortedParticlesI[--getParticleIndicesI[particles[i].cellIndex]] = particles[i];
    }
    copy(sortedParticlesI, sortedParticlesI + numParticles, particles);
}

void indexSortConstructionCountingSortImprovedParallel(Particle* particles, int numParticles, float h)
{
    boundingBoxConstruction(particles, numParticles, h);

    getParticleIndicesI.resize(boundingBox[6] + 1);
    memset(&getParticleIndicesI[0], 0, (boundingBox[6] + 1) * sizeof(getParticleIndicesI[0]));

    float invCellSize = 1.0f / (2.f * h);

    // #pragma omp parallel for
    for (int i = 0; i < numParticles; i++)
    {
        int k = static_cast<int>((particles[i].position.x - boundingBox[0]) * invCellSize);
        int l = static_cast<int>((particles[i].position.y - boundingBox[2]) * invCellSize);
        particles[i].k = k;
        particles[i].l = l;
        particles[i].cellIndex = k + l * boundingBox[4];
        getParticleIndicesI[particles[i].cellIndex] += 1;
    }

    for (size_t i = 1; i < boundingBox[6] + 1; i++)
    {
        getParticleIndicesI[i] += getParticleIndicesI[i - 1];
    }
    if (sortedParticlesI == nullptr)
    {
        sortedParticlesI = new Particle[numParticles];
    }
    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < numParticles; i++)
    {
        sortedParticlesI[--getParticleIndicesI[particles[i].cellIndex]] = particles[i];
    }
    copy(sortedParticlesI, sortedParticlesI + numParticles, particles);
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
                    for (size_t k = getParticleIndicesI.at(cellIndices[j]); k < getParticleIndicesI.at(cellIndices[j] + 1); k++)
                    {
                        Vector2f d = particles[i].position - particles[k].position;
                        float distance = d.x * d.x + d.y * d.y;
                        if (distance < (2.0f * h) * (2.0f * h))
                        {
                            particles[i].neighbors.push_back(k);
                        }
                    }
                }
            }
        }
    }
}

void indexSortQueryImproved(Particle* particles, int numParticles, float h) // too many adjacent cells at times (on the edges)
{
    // precompute
    float h2 = (2.0f * h) * (2.0f * h);

    // significant improvement by removing the offset calculation
    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            particles[i].neighbors.clear();
            // particles[i].neighbors.reserve(14); // only makes first iteration better
  
            int cellIndices[] = {
                particles[i].cellIndex,
                particles[i].cellIndex + 1,
                particles[i].cellIndex - 1,
                particles[i].cellIndex + boundingBox[4],
                particles[i].cellIndex - boundingBox[4],
                particles[i].cellIndex + boundingBox[4] + 1,
                particles[i].cellIndex + boundingBox[4] - 1,
                particles[i].cellIndex - boundingBox[4] + 1,
                particles[i].cellIndex - boundingBox[4] - 1
            };

            for (size_t j = 0; j < 9; j++)
            {
                if (cellIndices[j] >= 0 && cellIndices[j] < boundingBox[6])
                {
                    for (size_t k = getParticleIndicesI[cellIndices[j]]; k < getParticleIndicesI[cellIndices[j] + 1]; k++)
                    {
                        // this somehow is good improvement
                        float dx = particles[i].position.x - particles[k].position.x;
                        if (dx * dx >= h2) continue;

                        float dy = particles[i].position.y - particles[k].position.y;
                        float distance = dx * dx + dy * dy;
                        if (distance < h2) {
                            particles[i].neighbors.push_back(k);
                        }
                    }
                }
            }
        }
    }
}

void indexSortQueryKLImproved(Particle* particles, int numParticles, float h) // slightly worse
{
    // precompute
    float h2 = (2.0f * h) * (2.0f * h);
    int newIndexX;
    int newIndexY;

    // significant improvement by removing the offset calculation

    for (size_t i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            particles[i].neighbors.clear();
            // particles[i].neighbors.reserve(14); // only makes first iteration better

            for (size_t j = 0; j < 9; j++)
            {
                newIndexX = particles[i].k + cellOffset[j][0]; // preallocated cellOffset
                newIndexY = particles[i].l + cellOffset[j][1];

                if (newIndexX >= 0 && newIndexY >= 0 && newIndexX < boundingBox[4] && newIndexY < boundingBox[5])
                {
                    int index = newIndexX + newIndexY * boundingBox[4];
                    for (size_t k = getParticleIndicesI[index]; k < getParticleIndicesI[index + 1]; k++)
                    {
                        // this somehow is good improvement
                        float dx = particles[i].position.x - particles[k].position.x;
                        if (dx * dx >= h2) continue;

                        float dy = particles[i].position.y - particles[k].position.y;
                        float distance = dx * dx + dy * dy;
                        if (distance < h2) {
                            particles[i].neighbors.push_back(k);
                        }
                    }
                }
            }
        }
    }
}

void indexSortQueryOverCells(Particle* particles, int numParticles, float h)
{
    for (size_t i = 0; i < boundingBox[6]; i++)
    {
        // saving in terms of cell index computation does improve the performance slightly
        int cellIndices[] = {
            i,
            i + 1,
            i - 1,
            i + boundingBox[4],
            i - boundingBox[4],
            i + boundingBox[4] + 1,
            i + boundingBox[4] - 1,
            i - boundingBox[4] + 1,
            i - boundingBox[4] - 1
        };

        for (size_t j = getParticleIndicesI.at(i); j < getParticleIndicesI.at(i + 1); j++)
        {
            if (particles[j].isFluid)
            {
                particles[j].neighbors.clear();

                for (size_t k = 0; k < 9; k++)
                {
                    if (cellIndices[k] >= 0 && cellIndices[k] < boundingBox[4] * boundingBox[5])
                    {
                        for (size_t y = getParticleIndicesI.at(cellIndices[k]); y < getParticleIndicesI.at(cellIndices[k] + 1); y++)
                        {
                            Vector2f d = particles[j].position - particles[y].position;
                            float distance = d.x * d.x + d.y * d.y;
                            if (distance < (2.0f * h) * (2.0f * h))
                            {
                                particles[j].neighbors.push_back(y);
                            }
                        }
                    }
                }
            }
        }
    }
}

void indexSortQueryOverCellsImproved(Particle* particles, int numParticles, float h)
{
    // slow but probably nessecary
    float h2 = (2.0f * h) * (2.0f * h);

    for (size_t i = 0; i < boundingBox[6]; i++)
    {
        // saving in terms of cell index computation does improve the performance slightly
        // pretty similar speed in the improved versions
        int cellIndices[] = {
                    i,
                    i + 1,
                    i - 1,
                    i + boundingBox[4],
                    i - boundingBox[4],
                    i + boundingBox[4] + 1,
                    i + boundingBox[4] - 1,
                    i - boundingBox[4] + 1,
                    i - boundingBox[4] - 1
        };

        for (size_t j = getParticleIndicesI[i]; j < getParticleIndicesI[i + 1]; j++)
        {
            if (particles[j].isFluid)
            {
                particles[j].neighbors.clear();
                // particles[j].neighbors.reserve(14);

                for (size_t k = 0; k < 9; k++)
                {
                    if (cellIndices[k] >= 0 && i + cellIndices[k] < boundingBox[6])
                    {
                        for (size_t y = getParticleIndicesI[cellIndices[k]]; y < getParticleIndicesI[cellIndices[k] + 1]; y++)
                        {
                            float dx = particles[j].position.x - particles[y].position.x;
                            if (dx * dx >= h2) continue;

                            float dy = particles[j].position.y - particles[y].position.y;
                            float distance = dx * dx + dy * dy;
                            if (distance < h2) {
                                particles[j].neighbors.push_back(y);
                            }
                        }
                    }
                }
            }
        }
    }
}

void indexSortQueryCountNeighbors(Particle* particles, int numParticles, float h) // only count do not store
{
    if (numNeighbors == nullptr)
    {
        numNeighbors = new int[numParticles];
    }
    memset(&numNeighbors[0], 0, numParticles * sizeof(numNeighbors[0]));

    float h2 = (2.0f * h) * (2.0f * h);

    #pragma omp parallel for
    for (int i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            int cellIndices[] = {
                particles[i].cellIndex,
                particles[i].cellIndex + 1,
                particles[i].cellIndex - 1,
                particles[i].cellIndex + boundingBox[4],
                particles[i].cellIndex - boundingBox[4],
                particles[i].cellIndex + boundingBox[4] + 1,
                particles[i].cellIndex + boundingBox[4] - 1,
                particles[i].cellIndex - boundingBox[4] + 1,
                particles[i].cellIndex - boundingBox[4] - 1
            };

            for (size_t j = 0; j < 9; j++)
            {
                if (cellIndices[j] >= 0 && cellIndices[j] < boundingBox[6])
                {
                    for (size_t k = getParticleIndicesI[cellIndices[j]]; k < getParticleIndicesI[cellIndices[j] + 1]; k++)
                    {
                        float dx = particles[i].position.x - particles[k].position.x;
                        if (dx * dx >= h2) continue;

                        float dy = particles[i].position.y - particles[k].position.y;
                        float distance = dx * dx + dy * dy;
                        if (distance < h2) {
                            numNeighbors[i] += 1;
                        }
                    }
                }
            }
        }
    }
}

void indexSortQueryImprovedParallel(Particle* particles, int numParticles, float h)
{
    indexSortQueryCountNeighbors(particles, numParticles, h);

    float h2 = (2.0f * h) * (2.0f * h);

    #pragma omp parallel for
    for (int i = 0; i < numParticles; i++)
    {
        if (particles[i].isFluid)
        {
            int counter = numNeighbors[i];
            particles[i].neighbors.clear();
            particles[i].neighbors.resize(counter);
            counter--;

            int cellIndices[] = {
                particles[i].cellIndex,
                particles[i].cellIndex + 1,
                particles[i].cellIndex - 1,
                particles[i].cellIndex + boundingBox[4],
                particles[i].cellIndex - boundingBox[4],
                particles[i].cellIndex + boundingBox[4] + 1,
                particles[i].cellIndex + boundingBox[4] - 1,
                particles[i].cellIndex - boundingBox[4] + 1,
                particles[i].cellIndex - boundingBox[4] - 1
            };

            for (size_t j = 0; j < 9; j++)
            {
                if (cellIndices[j] >= 0 && cellIndices[j] < boundingBox[6])
                {
                    for (size_t k = getParticleIndicesI[cellIndices[j]]; k < getParticleIndicesI[cellIndices[j] + 1]; k++)
                    {
                        float dx = particles[i].position.x - particles[k].position.x;
                        if (dx * dx >= h2) continue;

                        float dy = particles[i].position.y - particles[k].position.y;
                        float distance = dx * dx + dy * dy;
                        if (distance < h2) {
                            particles[i].neighbors[counter] = k;
                            counter--;
                        }
                    }
                }
            }
        }
    }
}