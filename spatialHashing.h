#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"
#include "helperFunctions.h"

using namespace std;

/* */
void spatialHashingConstruction(Particle* particles, int numParticles, float h);
/* */
void spatialHashingQuery(Particle* particles, int numFluidParticles, float h);