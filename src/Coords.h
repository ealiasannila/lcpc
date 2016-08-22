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
#include <map>
class Coords {
private:
	std::map<int,std::set<Coords*>> leftNeighbours;
	std::map<int,std::set<Coords*>> rightNeighbours;

	double x;
	double y;

public:
	Coords();
	Coords(double x, double y);
	Coords(double x, double y, int polygon);
	double getX(){return x;}
	double getY(){return y;}
	int isRight(Coords* c1, Coords* c2);
	double eucDist(Coords* c1);
	double eucDistSquared(Coords* c1);
	std::pair<std::set<Coords*>::iterator, std::set<Coords*>::iterator> getRightNeighbours(int polygon);
	std::pair<std::set<Coords*>::iterator, std::set<Coords*>::iterator> getLeftNeighbours(int polygon);
	std::pair<std::map<int,std::set<Coords*>>::iterator,std::map<int,std::set<Coords*>>::iterator> belongsToPolygons();
	void addToPolygon(int polygon);
	void addNeighbours(Coords* l, Coords* r, int polygon);
	std::string toString();
	virtual ~Coords();
};

#endif /* COORDS_H_ */
