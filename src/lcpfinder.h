/*
 * lcpfinder.h
 *
 *  Created on: Aug 25, 2016
 *      Author: elias
 */

#ifndef SRC_LCPFINDER_H_
#define SRC_LCPFINDER_H_
#include <algorithm>
#include <iostream>
#include "Coords.h"
#include "Funnel.h"
#include "../lib/geometry.h"
#include <tr1/functional>
#include <tr1/unordered_set>

struct CoordsHasher {
	std::size_t operator()(const Coords& c) const {
		std::size_t res = 17;
		res = res * 31 + std::tr1::hash<double>()(c.getX());
		res = res * 31 + std::tr1::hash<double>()(c.getY());
		return res;
	}
};
class LcpFinder{
private:
	std::tr1::unordered_set<Coords, CoordsHasher> coordmap;
	std::map<int, double> frictions; //update to something ligther
	const Coords* getOpposing(const Coords* l, const Coords* r, int polygon);
	std::deque<Funnel> initFQue(const Coords* c, int polygon, nSet*neighbours);
	void findNeighboursInPolygon(const Coords* c, int polygon, nSet* neighbours);
	nSet findNeighbours(const Coords* c);
public:
	std::vector<Coords> leastCostPath(const Coords s, const Coords e);
	void addPolygon(int polygon, std::vector<std::vector<Coords>> points, double friction);

};



#endif /* SRC_LCPFINDER_H_ */
