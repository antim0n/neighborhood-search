#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include "particle.h"
#include "helperFunctions.h"

using namespace std;

/* */
void gridConstruction(Particle* particles, int numParticles, float h);
/* */
void gridConstructionGenerallyImproved(Particle* particles, int numParticles, float h);
/* */
void gridQuery(Particle* particles, int numFluidParticles, float h);
/* */
void gridQueryOverCells(float h);
/* */
void gridQueryGenerallyImproved(Particle* particles, int numFluidParticles, float h);