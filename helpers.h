#ifndef HELPERS_H
#define HELPERS_H

#include <random>
#include "vector2d.h"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> random_generator(0.,1.);

Vec2D random2D(Vec2D minDom, Vec2D maxDom)
{
    Vec2D deltaDom = maxDom - minDom;
    return Vec2D(random_generator(gen)*(deltaDom.X) + minDom.X , random_generator(gen)*(deltaDom.Y) + minDom.Y);
};

double random_local(double a, double b)
{
    return random_generator(gen)*(b - a) + a;
};

#endif