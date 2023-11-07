#ifndef IO_H
#define IO_H

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "Particle.h"
#include "wall.h"
#include "vector2d.h"


void readPositions(std::string nameFile, std::vector<Particle> *particles);

void readFixedParticles(std::string nameFile, std::vector<Particle> *particles);

void readRadii(std::string nameFile, std::vector<Particle> *particles);

void readMasses(std::string nameFile, std::vector<Particle> *particles);

void readWalls(std::string nameFile, std::vector<Wall> *walls);


#endif