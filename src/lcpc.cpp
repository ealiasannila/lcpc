//============================================================================
// Name        : lcpc.cpp
// Author      : Elias Annila
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "Coords.h"
#include "Funnel.h"
#include "minHeap.h"
#include "../lib/geometry.h"
#include <algorithm>
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

void printN(nContainer n) {
	for (const Coords* c : n) {
		std::cout << c->toString() << std::endl;
	}
}

/*
 * Finds opposing by looking up intersection of the immidiate neighbours of either end of the base.
 * This can be either 1 or 2 vertices, of whcih 1 or 0 may be relevant. (irrelevant vertice is already inside funnel)
 * LOOKING UP INTERSECTION REQUIRES SORTING IMMIDIATE NEIGHBOURS. HOW TO MAKE FASTER:
 * ?-> have sorted vector of immidiate neighbours available. (increases insertion time to log N, now constant)
 */
const Coords* getOpposing(const Coords* l, const Coords* r, int polygon) {
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

std::deque<Funnel> initFQue(const Coords* c, int polygon, std::set<const Coords*>*neighbours) {
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

void findNeighboursInPolygon(const Coords* c, int polygon, std::set<const Coords*>* neighbours) {
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
std::set<const Coords*> findNeighbours(const Coords* c) {
	std::set<const Coords*> neighbours { };
	allNeighIter polyIt = c->getAllLeftN(); //Only keys area interesting here, first = begin iterator, second end iterator
	for (allNContainer::iterator it = polyIt.first; it != polyIt.second; it++) {
		findNeighboursInPolygon(c, it->first, &neighbours);
	}
	return neighbours;
}

void leastCostPath(const Coords* start, const Coords* end) {
	start->setToStart(0);
	struct coordComp {
			bool operator()(const Coords* x, const Coords* y) const {
				return ((x->getToStart())>(y->getToStart()));
			}
	} compare { };
	MinHeap<Coords*, coordComp> minheap(compare);

	std::vector<Coords> cv { Coords(1,1),Coords(1,2),Coords(1,3),Coords(1,4) } ;
	for (int i = 0; i<cv.size(); i++ ) {
		cv[i].setToStart(i);
		minheap.push(&cv[i]);
	}
	while (!minheap.empty()) {
		std::cout << minheap.top()->toString()<<"tostart: "<<minheap.top()->getToStart() <<std::endl;
		minheap.pop();
	}
	//MinHeap<const Coords*> minHeap;
	//minHeap.push(start);

}

int main() {

	Polygon poly("sample1.bdm", true);
	poly.setDebugOption(false);      //set debug flag;
	poly.triangulation();
	std::cout << "triangulaton done" << std::endl;
	PointbaseMap points = poly.points();
	std::tr1::unordered_set<Coords, CoordsHasher> coordmap;

	Triangles triangles = poly.triangles();

	int t = 1;
	for (std::list<Triangle>::iterator it = triangles.begin(); it != triangles.end(); it++) {
		std::cout << "\n\n--Triangle: " << t << "--\n";
		t++;
		Triangle triangle = *it;
		const Coords* cp[3];
		Coords c[3];
		for (unsigned i = 0; i < 3; i++) {
			Pointbase pb = *points.at(triangle[i]);
			c[i] = Coords(pb.x, pb.y, 0);

		}
		// if triangle orientation is clockwise turn it to CCW
		if (c[0].isRight(&c[1], &c[2]) == 1) {
			Coords tmp = c[2];
			c[2] = c[1];
			c[1] = tmp;
			std::cout << "changed orientation." << std::endl;
		} else {
			std::cout << "orientation fine" << std::endl;
		}

		for (unsigned i = 0; i < 3; i++) {
			std::pair<std::tr1::unordered_set<Coords>::iterator, bool> p = coordmap.insert(c[i]);
			cp[i] = &*p.first;
			bool alreadyThere = !p.second;
			if (alreadyThere) {
				std::cout << "merging at " << c[i].toString() << std::endl;
			} else {
				std::cout << "inserted new at " << c[i].toString() << std::endl;
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
			std::cout << "adding to" << cp[i]->toString() << std::endl;
			//std::cout<<"L"<<->toString()<<std::endl;
			//std::cout<<"R"<<&*coordmap.find(c[r])->toString()<<std::endl;
			cp[i]->addNeighbours(cp[l], cp[r], 0);
		}
	}

	std::cout << "coords found: " << coordmap.size() << std::endl;
//set start vertex:
	std::tr1::unordered_set<Coords>::iterator sit = coordmap.find(Coords(400, 50, 0));
	const Coords* s = &*sit;
	std::tr1::unordered_set<Coords>::iterator fit = coordmap.find(Coords(50, 200, 0));
	const Coords* f = &*fit;

	leastCostPath(s, f);
	return 0;
}

