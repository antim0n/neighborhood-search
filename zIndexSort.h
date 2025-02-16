#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"
#include "helperFunctions.h"

using namespace std;

// z-index sort
extern Handle* sortedIndicesZI;

/* */
void zIndexSortConstruction(Particle* particles, int numParticles, float h);
/* */
void zIndexSortConstructionImproved(Particle* particles, int numParticles, float h);
/* */
void zIndexSortConstructionHandleSort(Particle* particles, int numParticles, float h);
/* */
void zIndexSortConstructionHandleSortImproved(Particle* particles, int numParticles, float h);
/* */
void zIndexSortQuery(Particle* particles, int numParticles, float h);
/* */
void zIndexSortQueryImproved(Particle* particles, int numParticles, float h);
/* */
void zIndexSortQueryHandleSort(Particle* particles, int numParticles, float h);
/* */
void zIndexSortQueryHandleSortImproved(Particle* particles, int numParticles, float h);
/* */
void zIndexSortQueryHandleSortOverCellsImproved(Particle* particles, int numParticles, float h);