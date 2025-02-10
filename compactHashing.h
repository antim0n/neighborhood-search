#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"
#include "helperFunctions.h"

using namespace std;

extern Handle* sortedIndicesCH;

// compact hashing
// extern int* handleArray;
// extern vector<Particle*>* compactList;
// extern vector<Particle*>* usedCellsEntries;

/* */
void compactHashingConstruction(Particle* particles, int numParticles, float h);
/* */
void compactHashingConstructionZSorted(Particle* particles, int numParticles, float h);
/* */
void compactHashingQuery(Particle* particles, int numParticles, float h);