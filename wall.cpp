#include "wall.h"

Wall::Wall(){
	wall.start.set(0.0, 0.0);
	wall.end.set(0.0, 0.0);
}

Wall::Wall(double x1, double y1, double x2, double y2) {
	wall.start.set(x1, y1);
	wall.end.set(x2, y2);
}

Wall::Wall(Vec2D start, Vec2D end) {
	wall.start = start;
	wall.end = end;
}


Vec2D Wall::getNearestPoint(Vec2D position_i) {
	Vec2D relativeEnd, relativePos, relativeEndScal, relativePosScal;
	double dotProduct;
	Vec2D nearestPoint;

	// Create Vector Relative to Wall's 'start'
	relativeEnd = wall.end - wall.start;	// Vector from wall's 'start' to 'end'
	relativePos = position_i - wall.start;	// Vector from wall's 'start' to agent i 'position'

	// Scale Both Vectors by the Length of the Wall
	relativeEndScal = relativeEnd;
	relativeEndScal.normalise();

	relativePosScal = relativePos * (1.0F / relativeEnd.getLength());

	// Compute Dot Product of Scaled Vectors
	dotProduct = Vec2D::dot(relativeEndScal,relativePosScal);

	if (dotProduct < 0.0)		// Position of Agent i located before wall's 'start'
    {
		nearestPoint = wall.start;
    }
    else if (dotProduct > 1.0)	// Position of Agent i located after wall's 'end'
		nearestPoint = wall.end;
	else						// Position of Agent i located between wall's 'start' and 'end'
		nearestPoint = (relativeEnd * dotProduct) + wall.start;

	return nearestPoint;
}

Vec2D Wall::getNormal(){
    return normal;
}

int Wall::getIndex()
{
    return index;
}

void Wall::setNormal(Vec2D normalIn){
    normal = Vec2D::getUnitVector(normalIn);
}

void Wall::setIndex(int indexWall){
    index = indexWall;
}
   
