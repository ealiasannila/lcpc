/*
 * Coords.h
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#ifndef COORDS_H_
#define COORDS_H_

class Coords {
private:
	double x;
	double y;
public:
	Coords(double x, double y);
	double getX(){return x;}
	double getY(){return y;}
	virtual ~Coords();
};

#endif /* COORDS_H_ */
