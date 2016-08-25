/*
 * lcpfinder.cpp
 *
 *  Created on: Aug 25, 2016
 *      Author: elias
 */


#include "lcpfinder.h"
#include "minHeap.h"

/*
 * Finds opposing by looking up intersection of the immidiate neighbours of either end of the base.
 * This can be either 1 or 2 vertices, of whcih 1 or 0 may be relevant. (irrelevant vertice is already inside funnel)
 * LOOKING UP INTERSECTION REQUIRES SORTING IMMIDIATE NEIGHBOURS. HOW TO MAKE FASTER:
 * ?-> have sorted vector of immidiate neighbours available. (increases insertion time to log N, now constant)
 */
const Coords* LcpFinder::getOpposing(const Coords* l, const Coords* r, int polygon) {
	nContainer ln = l->getLeftNeighbours(polygon);
	nContainer rn = r->getRightNeighbours(polygon);

	sort(ln.begin(), ln.end()); //... could this be avoided? Or is it that bad? they are small vectors...
	sort(rn.begin(), rn.end()); //...

	nContainer intersection;
	std::set_intersection(ln.begin(), ln.end(), rn.begin(), rn.end(), std::back_inserter(intersection));
	for (nContainer::iterator it = intersection.begin(); it != intersection.end(); it++) {
		if (*it != 0 and (*it)->isRight(l, r) == -1) {
			return *it;
		}
	}
	return 0;
}

std::deque<Funnel> LcpFinder::initFQue(const Coords* c, int polygon, std::set<const Coords*>*neighbours) {
	nContainer ln = c->getLeftNeighbours(polygon);
	nContainer rn = c->getRightNeighbours(polygon);

	std::deque<Funnel> funnelQue { };
	nContainer::iterator rit = rn.begin();
	for (nContainer::iterator lit = ln.begin(); lit != ln.end(); lit++) {
		Funnel f(*lit, c, *rit);
		neighbours->insert(*lit);
		neighbours->insert(*rit);
		funnelQue.push_back(f);
		rit++;
	}
	return funnelQue;

}

void LcpFinder::findNeighboursInPolygon(const Coords* c, int polygon, nSet* neighbours) {
	std::deque<Funnel> funnelQue = initFQue(c, polygon, neighbours);

	while (!funnelQue.empty()) {
		Funnel f = funnelQue.front();
		funnelQue.pop_front();
		std::pair<const Coords*, const Coords*> base = f.getBase();
		const Coords* o = getOpposing(base.first, base.second, polygon);
		if (o != 0) {
			f.reactToOpposite(o, &funnelQue, neighbours);
			funnelQue.push_back(f);
		}
	}
}
nSet LcpFinder::findNeighbours(const Coords* c) {
	nSet neighbours { };
	allNeighIter polyIt = c->getAllLeftN(); //Only keys area interesting here, first = begin iterator, second end iterator
	for (allNContainer::iterator it = polyIt.first; it != polyIt.second; it++) {
		findNeighboursInPolygon(c, it->first, &neighbours);
	}
	return neighbours;
}

void LcpFinder::leastCostPath(const Coords* start, const Coords* end) {
	start->setToStart(0);
	struct compToStart {
		bool operator()(const Coords* x, const Coords* y) const {
			return ((x->getToStart()) > (y->getToStart()));
		}
	} compare { };
	MinHeap<const Coords*, compToStart> minheap(compare);

	minheap.push(start);

	while (!minheap.empty()) {
		const Coords* node = minheap.top();
		minheap.pop();
		nSet neighbours = findNeighbours(node);
		for (const Coords* n : neighbours) {
			double d { node->getToStart() + node->eucDist(n) };
			if (n->getToStart() < 0) { // node has not yet been inserted into minheap
				n->setToStart(d);
				n->setPred(node);
				minheap.push(n);
			} else if (n->getToStart() > d) {
				n->setToStart(d);
				n->setPred(node);
				minheap.update(); //reorders minheap after changing priority of n
			}
		}
	}
	//MinHeap<const Coords*> minHeap;
	//minHeap.push(start);
	std::cout << "ROUTE FROM END TO START:" << std::endl;
	const Coords* pred = end->getPred();
	std::cout << end->toString() << std::endl;
	while (pred != 0) {
		std::cout << pred->toString() << std::endl;
		pred = pred->getPred();
	}

}

void LcpFinder::addPolygon(int polygon, std::vector<std::vector<Coords>> points) {
	for (std::vector<Coords> ring : points) {
		for (Coords p : ring) {
			p.addToPolygon(polygon);
			// SOMEHOW ADD TO TRIANGULATOR...
		}
	}
	/*
	 * for triangle in triangles:
	 * 		this->coordmap.insert(p);
	 *
	 */

}

