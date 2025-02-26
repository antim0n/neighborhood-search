#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"
#include "helperFunctions.h"

using namespace std;

extern int* compactCellIndex;
extern int* compactReference;

/* */
void cellLinkedListConstruction(Particle* particles, int numParticles, float h);

/* */
void cellLinkedListQuery(Particle* particles, int numParticles, float h);