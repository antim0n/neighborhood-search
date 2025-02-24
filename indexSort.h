#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"
#include "helperFunctions.h"

using namespace std;

extern Particle* sortedParticlesI;

/* */
void indexSortConstructionCountingSort(Particle* particles, int numParticles, float h);
/* */
void indexSortConstructionCompareSorting(Particle* particles, int numParticles, float h);
/* */
void indexSortConstructionCountingSortImproved(Particle* particles, int numParticles, float h);
/* */
void indexSortConstructionInsertionSortImproved(Particle* particles, int numParticles, float h);
/* */
void indexSortConstructionCountingSortImprovedParallel(Particle* particles, int numParticles, float h);

/* */
void indexSortQuery(Particle* particles, int numParticles, float h);
/* */
void indexSortQueryImproved(Particle* particles, int numParticles, float h);
/* */
void indexSortQueryKLImproved(Particle* particles, int numParticles, float h);
/* */
void indexSortQueryOverCells(Particle* particles, int numParticles, float h);
/* */
void indexSortQueryOverCellsImproved(Particle* particles, int numParticles, float h);
/* */
void indexSortQueryImprovedParallel(Particle* particles, int numParticles, float h);