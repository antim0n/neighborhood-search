#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"

using namespace std;

// index sort
extern float* boundingBox; // xmin, xmax, ymin, ymax, cellsx, cellsy
extern vector<int> getParticleIndices; // change to vector, no memory allocations

// spatial hashing

/* */
void gridConstruction();
/* */
void gridQuery();
/* */
void indexSortConstruction(Particle* particles, int numFluidParticles, float H);
/* */
void indexSortQuery(Particle* particles, int numFluidParticles, float H);
/* */
void zIndexSortConstruction();
/* */
void zIndexSortQuery();
/* */
void spatialHashingConstruction();
/* */
void spatialHashingQuery();
/* */
void compactHashingConstruction();
/* */
void compactHashingQuery();