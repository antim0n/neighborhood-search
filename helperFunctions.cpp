#include "helperFunctions.h"

float* boundingBox = new float[6];

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

unsigned int hashFunction(int cellIndexX, int cellIndexY, unsigned int sizeHashTable)
{
    return ((cellIndexX * 73856093u) ^ (cellIndexY * 19349663u)) % sizeHashTable; // overflow creates negative number, either add sizeHashTable, use unsigned int or other larger datatypes
}

// Helper function to spread the bits of a 32-bit integer into 64 bits // https://lemire.me/blog/2018/01/08/how-fast-can-you-bit-interleave-32-bit-integers/
uint64_t spreadBits(uint32_t x) {
    uint64_t result = x;
    result = (result | (result << 16)) & 0x0000FFFF0000FFFF;
    result = (result | (result << 8)) & 0x00FF00FF00FF00FF;
    result = (result | (result << 4)) & 0x0F0F0F0F0F0F0F0F;
    result = (result | (result << 2)) & 0x3333333333333333;
    result = (result | (result << 1)) & 0x5555555555555555;
    return result;
}

// Function to interleave bits of two 32-bit integers x and y
uint64_t interleaveBits(uint32_t x, uint32_t y) {
    return (spreadBits(x) | (spreadBits(y) << 1));
}