/*
 * const Coords.h
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#ifndef  Coords_H_
#define  Coords_H_
#include <string>
#include <set>
#include <map>
#include <cmath>
class  Coords {
private:
	mutable std::map<int, std::set< const Coords*>> leftNeighbours;
	mutable std::map<int, std::set< const Coords*>> rightNeighbours;

	double x;
	double y;

public:
	 Coords();
	 Coords(double x, double y);
	 Coords(double x, double y, int polygon);
	double getX()  const{
		return x;
	}
	double getY()  const{
		return y;
	}
	int isRight(const Coords* c1, const Coords* c2) const;
	double eucDist(const Coords* c1);
	double eucDistSquared(const Coords* c1);
	std::pair<std::set< const Coords*>::iterator, std::set< const Coords*>::iterator> getRightNeighbours(int polygon) const;
	std::pair<std::set< const Coords*>::iterator, std::set< const Coords*>::iterator> getLeftNeighbours(int polygon) const;
	std::pair<std::map<int, std::set< const Coords*>>::iterator, std::map<int, std::set< const Coords*>>::iterator> getAllLeftN() const;
	std::pair<std::map<int, std::set< const Coords*>>::iterator, std::map<int, std::set< const Coords*>>::iterator> getAllRightN() const;
	void addToPolygon(int polygon)const ;
	void addNeighbours(const Coords* l,const Coords* r, int polygon) const;
	std::string toString() const;
	virtual ~ Coords();

	bool operator==( const Coords& c) const;
	bool operator!=( const Coords& c) const;
};

#endif /* const Coords_H_ */
