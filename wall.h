#ifndef WALL_H
#define WALL_H

#include "vector2d.h"
#include "Particle.h"

struct Line {
	Vec2D start;
	Vec2D end;
};

class Wall {
private:
	Line wall;
    Vec2D normal;
    int index;

public:
	Wall();
	Wall(double x1, double y1, double x2, double y2);
    Wall(Vec2D start, Vec2D end);
	//~Wall();

	Vec2D getStartPoint() const { return wall.start; }
	Vec2D getEndPoint() const { return wall.end; }
	Vec2D getNearestPoint(Vec2D position_i);	// Computes distance between 'position_i' and wall
    Vec2D getNormal();
    void setIndex(int indexWall);
    int getIndex();
    void setNormal(Vec2D normalIn);
};

#endif