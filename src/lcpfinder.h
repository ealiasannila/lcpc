/*
 * lcpfinder.h
 *
 *  Created on: Aug 25, 2016
 *      Author: elias
 */

#ifndef SRC_LCPFINDER_H_
#define SRC_LCPFINDER_H_
#include <unordered_set>
struct CoordsHasher;
class LcpFinder{
private:
	std::tr1::unordered_set<Coords, CoordsHasher> coordmap;

	const Coords* getOpposing(const Coords* l, const Coords* r, int polygon);
	std::deque<Funnel> initFQue(const Coords* c, int polygon, std::set<const Coords*>*neighbours);
	void findNeighboursInPolygon(const Coords* c, int polygon, nSet* neighbours);
	nSet findNeighbours(const Coords* c);
public:
	void leastCostPath(const Coords* start, const Coords* end);
	void addPolygon(int polygon, std::vector<std::vector<Coords>> points);

};



#endif /* SRC_LCPFINDER_H_ */
