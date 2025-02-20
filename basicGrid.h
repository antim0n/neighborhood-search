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
void gridConstructionImproved(Particle* particles, int numParticles, float h);
/* */
void gridConstructionImprovedParallel(Particle* particles, int numParticles, float h);
/* */
void gridQuery(Particle* particles, int numFluidParticles, float h);
/* */
void gridQueryOverCells(Particle* particles, float h);
/* */
void gridQueryImproved(Particle* particles, int numFluidParticles, float h);
/* */
void gridQueryImprovedParallel(Particle* particles, int numFluidParticles, float h);