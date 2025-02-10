#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"
#include "helperFunctions.h"

using namespace std;

/* */
void indexSortConstruction(Particle* particles, int numParticles, float h);
/* */
void indexSortQuery(Particle* particles, int numParticles, float h);