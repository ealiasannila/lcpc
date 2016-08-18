/*
 * Coords.h
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#ifndef COORDS_H_
#define COORDS_H_
#include <string>

class Coords {
private:
	double x;
	double y;
public:
	Coords(double x, double y);
	double getX(){return x;}
	double getY(){return y;}
	int isRight(Coords* c1, Coords* c2);
	double eucDist(Coords* c1);
	double eucDistSquared(Coords* c1);
	std::string toString();
	virtual ~Coords();
};

#endif /* COORDS_H_ */
