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

std::vector<Coords> LcpFinder::leastCostPath(const Coords s, const Coords e) {

	std::tr1::unordered_set<Coords>::iterator startIt = this->coordmap.find(s);
	const Coords* start;
	if (startIt != coordmap.end()) {
		start = &*startIt;
	}
	std::tr1::unordered_set<Coords>::iterator endIt = this->coordmap.find(e);
	const Coords* end;
	if (endIt != coordmap.end()) {
		end = &*endIt;
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
	std::cout<<"nearly\n";
	std::vector<Coords> res { *end };
	std::cout<<"almost\n";
	return res;
}

std::vector<Coords> LcpFinder::intermidiatePoints(Coords p, Coords next, double maxDist) {
	std::cout << "P: " << p.toString() << std::endl;
	std::cout << "N: " << next.toString() << std::endl;
	int pointsToAdd = std::ceil(p.eucDist(&next) / maxDist) - 1;
	std::cout << "adding: " << pointsToAdd << std::endl;
	std::vector<Coords> res(pointsToAdd);
	for (int i = 0; i < pointsToAdd; i++) {
		double x { ((next.getX() - p.getX()) / (pointsToAdd + 1)) * (i + 1) + p.getX() };
		double y { ((next.getY() - p.getY()) / (pointsToAdd + 1)) * (i + 1) + p.getY() };
		res[i] = Coords(x, y);
		std::cout << res[i].toString() << std::endl;
	}
	return res;
}

void LcpFinder::addPolygon(int polygon, std::vector<std::vector<Coords>> points, double friction, double maxDist) {
	this->frictions[polygon] = friction;
	unsigned int i { 1 };
	Polygon triangulator { };
	for (unsigned int ring = 0; ring < points.size(); ring++) {
		std::cout<<"ringsize: "<<points[ring].size()<<std::endl;
		for (unsigned int pi = 0; pi < points[ring].size(); pi++) {
			Coords p = points[ring][pi];
			p.addToPolygon(polygon);
			triangulator.addPoint(i, p.getX(), p.getY(), ring);
			i++;

			// add intermidiate points if next point is too far.
			Coords next = (pi != points[ring].size() - 1) ? points[ring][pi + 1] : points[ring][0];
			if (next.eucDist(&p)+0.0000001 > maxDist) {
				for (Coords ip : this->intermidiatePoints(p, next, maxDist)) {
					std::cout << "ip: " << ip.getX() << ":" << ip.getY() << std::endl;
					ip.addToPolygon(polygon);
					triangulator.addPoint(i, ip.getX(), ip.getY(), ring);
					i++;
				}
			}

		}
	}
	std::cout << "points added\n";
	triangulator.addEdges();
	std::cout << "edges added\n";
	std::cout << triangulator.points().at(1)->id;
	triangulator.initializate();
	std::cout << "init done\n";
	triangulator.triangulation();
	Triangles triangles = triangulator.triangles();
	triangulator.saveAsShowme();
	std::cout << "triangulaton done: " << triangles.size() << " triangles created" << std::endl;



	for (std::list<Triangle>::iterator it = triangles.begin(); it != triangles.end(); it++) {
		Triangle triangle = *it;
		const Coords* cp[3];
		Coords c[3];
		for (unsigned i = 0; i < 3; i++) {
			Pointbase pb = *triangulator.points().at(triangle[i]);
			c[i] = Coords(pb.x, pb.y, polygon);

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
				cp[i]->addToPolygon(polygon);
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
			cp[i]->addNeighbours(cp[l], cp[r], polygon);
		}
	}

}

