/*
 * lcpfinder.cpp
 *
 *  Created on: Aug 25, 2016
 *      Author: elias
 */

#include "lcpfinder.h"
#include "minHeap.h"
#include <array>

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

std::deque<Funnel> LcpFinder::initFQue(const Coords* c, int polygon, nSet*neighbours) {
	nContainer ln = c->getLeftNeighbours(polygon);
	nContainer rn = c->getRightNeighbours(polygon);

	std::deque<Funnel> funnelQue { };
	nContainer::iterator rit = rn.begin();
	for (nContainer::iterator lit = ln.begin(); lit != ln.end(); lit++) {
		Funnel f(*lit, c, *rit);
		neighbours->insert(std::pair<const Coords*, int>(*lit, polygon));
		neighbours->insert(std::pair<const Coords*, int>(*rit, polygon));
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
			f.reactToOpposite(o, &funnelQue, neighbours, polygon);
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

std::vector<Coords> LcpFinder::leastCostPath(Coords s, std::vector<Coords> e) {

	std::tr1::unordered_set<Coords>::iterator startIt = this->coordmap.find(s);
	const Coords* start;
	if (startIt != coordmap.end()) {
		start = &*startIt;
	}


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
		for (std::pair<const Coords*, int> p : neighbours) {
			const Coords* n = p.first;
			double d { node->getToStart() + node->eucDist(n) * this->frictions[p.second] };
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
	//update to many endpoints
	std::vector<Coords> res;
	for(Coords ep : e){
		std::tr1::unordered_set<Coords>::iterator endit = this->coordmap.find(ep);
		if(endit != coordmap.end()){
			res.push_back(*endit);
		}
	}
	return res;
}


void LcpFinder::addPolygon(int polygon, std::vector<p2t::Point*> steinerPoints, std::vector<std::vector<p2t::Point*>> points, double friction) {
	this->frictions[polygon] = friction;
		p2t::CDT cdt { points[0] }; //Constrained Delaunay triangulator with outer edges.

	for (unsigned int hole = 1; hole < points.size(); hole++) {
		cdt.AddHole(points[hole]);
	}
	for(p2t::Point* sp: steinerPoints){
		cdt.AddPoint(sp);
	}

	cdt.Triangulate();
	std::vector<p2t::Triangle*> triangles = cdt.GetTriangles();
	for (std::vector<p2t::Triangle*>::iterator it = triangles.begin(); it != triangles.end(); it++) {
		p2t::Triangle* triangle = *it;
		const Coords* cp[3];
		Coords c[3];
		for (unsigned i = 0; i < 3; i++) {
			p2t::Point* point = triangle->GetPoint(i);
			c[i] = Coords(point->x, point->y, polygon);

		}
		// if triangle orientation is clockwise turn it to CCW
		if (c[0].isRight(&c[1], &c[2]) == 1) {
			Coords tmp = c[2];
			c[2] = c[1];
			c[1] = tmp;
		}
		for (unsigned i = 0; i < 3; i++) {
			std::pair<std::tr1::unordered_set<Coords>::iterator, bool> p = coordmap.insert(c[i]);
			cp[i] = &*p.first;
			bool alreadyThere = !p.second;
			if (alreadyThere) {
				cp[i]->addToPolygon(polygon);
			}
		}
		for (unsigned i = 0; i < 3; i++) {
			int l = i - 1;
			int r = i + 1;
			if (l == -1) {
				l = 2;
			}
			if (r == 3) {
				r = 0;
			}
			cp[i]->addNeighbours(cp[l], cp[r], polygon);
		}
	}

}

