#pragma once
#include "particle.h"

struct Handle
{
	int cellIndex;
	Particle* reference;
};

extern float* boundingBox; // xmin, xmax, ymin, ymax, cellsx, cellsy, cells

static int cellOffset[][2] = {
    {-1, 1}, {0, 1}, {1, 1},
    {-1, 0}, {0, 0}, {1, 0},
    {-1, -1}, {0, -1}, {1, -1}
};

/* */
void boundingBoxConstruction(Particle* particles, int numParticles, float h);
/* */
unsigned int hashFunction(int cellIndexX, int cellIndexY, unsigned int sizeHashTable);
/* */
uint64_t spreadBits(uint32_t x);
/* */
uint64_t interleaveBits(uint32_t x, uint32_t y);