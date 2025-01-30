#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"

using namespace std;

// grid
extern float* boundingBox; // xmin, xmax, ymin, ymax, cellsx, cellsy
extern vector<vector<Particle*>> cells;

// index sort
extern vector<int> getParticleIndices; // change to vector, no memory allocations

// spatial hashing
extern vector<Particle*> hashTable[];

// z-index sort
// extern int** sortedIndices;

// compact hashing
// extern int* handleArray;
// extern vector<Particle*>* compactList;
// extern vector<Particle*>* usedCellsEntries;

/* */
void boundingBoxConstruction(Particle* particles, int numParticles, float h);
/* */
void gridConstruction(Particle* particles, int numParticles, float h);
/* */
void gridQuery(Particle* particles, int numParticles, float h);
/* */
void indexSortConstruction(Particle* particles, int numParticles, float h);
/* */
void indexSortQuery(Particle* particles, int numParticles, float h);
/* */
void zIndexSortConstruction(Particle* particles, int numParticles, float h);
/* */
void zIndexSortQuery(Particle* particles, int numParticles, float h);
/* */
void spatialHashingConstruction(Particle* particles, int numParticles, float h);
/* */
void spatialHashingQuery(Particle* particles, int numFluidParticles, float h);
/* */
void compactHashingConstruction(Particle* particles, int numParticles, float h);
/* */
void compactHashingQuery(Particle* particles, int numParticles, float h);