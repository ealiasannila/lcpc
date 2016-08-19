/*
 * Coords.h
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#ifndef COORDS_H_
#define COORDS_H_
#include <string>
#include <set>
class Coords {
private:
	std::set<Coords*> leftNeighbours;
	std::set<Coords*> rightNeighbours;

	double x;
	double y;
public:
	Coords();
	Coords(double x, double y);
	double getX(){return x;}
	double getY(){return y;}
	int isRight(Coords* c1, Coords* c2);
	double eucDist(Coords* c1);
	double eucDistSquared(Coords* c1);
	std::set<Coords*>::iterator getLeftNeighbours();
	std::set<Coords*>::iterator getRightNeighbours();
	std::set<Coords*>::iterator getLeftNeighboursEnd();
	std::set<Coords*>::iterator getRightNeighboursEnd();
	void addNeighbours(Coords* l, Coords* r);
	std::string toString();
	virtual ~Coords();
};

#endif /* COORDS_H_ */
