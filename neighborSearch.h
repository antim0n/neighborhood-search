#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"

using namespace std;

// grid
extern vector<vector<Particle*>> cells;

// index sort
extern float* boundingBox; // xmin, xmax, ymin, ymax, cellsx, cellsy
extern vector<int> getParticleIndices; // change to vector, no memory allocations

// spatial hashing
extern vector<Particle*> hashTable[];

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
void zIndexSortConstruction();
/* */
void zIndexSortQuery();
/* */
void spatialHashingConstruction(Particle* particles, int numParticles, float h);
/* */
void spatialHashingQuery(Particle* particles, int numParticles, float h);
/* */
void compactHashingConstruction();
/* */
void compactHashingQuery();