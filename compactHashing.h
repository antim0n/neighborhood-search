#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"
#include "helperFunctions.h"

using namespace std;

extern int* hashTableCH;
extern Handle* sortedIndicesCH;

// compact hashing
// extern int* handleArray;
// extern vector<Particle*>* compactList;
// extern vector<Particle*>* usedCellsEntries;

/* */
void compactHashingConstruction(Particle* particles, int numParticles, float h);
/* */
void compactHashingConstructionImproved(Particle* particles, int numParticles, float h);
/* */
void compactHashingConstructionHashCollisionFlagImproved(Particle* particles, int numParticles, float h);
/* */
void compactHashingConstructionZSortedImproved(Particle* particles, int numParticles, float h);
/* */
void compactHashingConstructionHandleSort(Particle* particles, int numParticles, float h);
/* */
void compactHashingConstructionHandleSortImproved(Particle* particles, int numParticles, float h);
/* */
void compactHashingConstructionHandleSortImprovedParallel(Particle* particles, int numParticles, float h);

/* */
void compactHashingQuery(Particle* particles, int numParticles, float h);
/* */
void compactHashingQueryImproved(Particle* particles, int numParticles, float h);
/* */
void compactHashingQueryHashCollisionFlagImproved(Particle* particles, int numParticles, float h);
/* */
void compactHashingQueryImprovedParallel(Particle* particles, int numParticles, float h);