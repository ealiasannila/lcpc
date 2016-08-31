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
#include "minHeap.h"
#include "geomfunc.h"
#include "../lib/poly2tri.h"
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
	std::vector<std::vector<std::vector<p2t::Point*>>> polygons;
	std::map<int, std::vector<p2t::Point*>> targetPoints;
	std::vector<double> frictions;

	const Coords* getOpposing(const Coords* l, const Coords* r, int polygon);
	std::deque<Funnel> initFQue(const Coords* c, int polygon, nSet*neighbours);
	void findNeighboursInPolygon(const Coords* c, int polygon, nSet* neighbours);
	nSet findNeighbours(const Coords* c);
	std::vector<p2t::Point*> intermidiatePoints(p2t::Point p, p2t::Point next, double maxDist);
public:
	std::tr1::unordered_set<Coords, CoordsHasher> getCoordmap(){return coordmap;}
	std::vector<Coords> leastCostPath(Coords s, std::vector<Coords> e);
	void addPolygon(std::vector<std::vector<p2t::Point*>> points, double friction);
	void addSteinerPoints(std::vector<p2t::Point*> steinerpoints, int polygon);
	void addSteinerPoint(p2t::Point* steinerpoint, int polygon);
	void triangulate(int polygon);
    ~LcpFinder();
};



#endif /* SRC_LCPFINDER_H_ */
