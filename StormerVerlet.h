#ifndef STORMERVERLET_H
#define STORMERVERLET_H

#include "vector2d.h"
#include "Particle.h"

void updateX(Particle *p, double delta_t);

void updateV(Particle *p, double delta_t);

#endif