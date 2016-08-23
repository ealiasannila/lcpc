//============================================================================
// Name        : lcpc.cpp
// Author      : Elias Annila
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "Coords.h"
#include "Funnel.h"
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

const Coords* getOpposing(const Coords* l, const Coords* r, int polygon) {
	std::pair<std::set<const Coords*>::iterator, std::set<const Coords*>::iterator> ln = l->getLeftNeighbours(polygon);
	std::pair<std::set<const Coords*>::iterator, std::set<const Coords*>::iterator> rn = r->getRightNeighbours(polygon);
	/*std::cout<<"LEFT NEIGHBOURS OF 50-200\n";
	 for (std::set<const Coords*>::iterator it = ln.first; it != ln.second; it++) {
	 std::cout << (*it)->toString() << std::endl;
	 }
	 std::cout<<"RIGHT NEIGHBOURS OF 115-160\n";
	 for (std::set<const Coords*>::iterator it = rn.first; it != rn.second; it++) {
	 std::cout << (*it)->toString() << std::endl;
	 }
	 */
	std::vector<const Coords*> intersection;
	std::set_intersection(ln.first, ln.second, rn.first, rn.second, std::back_inserter(intersection));
	if (intersection.size() > 0) {
		return *intersection.begin();
	}

	std::cout << "No opposing found" << std::endl;
	return 0;
}

void findNeighboursInPolygon(const Coords* c, int polygon, std::list<const Coords*>* neighbours) {
	std::pair<std::set<const Coords*>::iterator, std::set<const Coords*>::iterator> ln = c->getLeftNeighbours(polygon);
	std::pair<std::set<const Coords*>::iterator, std::set<const Coords*>::iterator> rn = c->getRightNeighbours(polygon);
	std::cout << std::endl << std::endl;
	std::deque<Funnel> funnelQue { };
	std::set<const Coords*>::iterator rit = rn.first;


	for (std::set<const Coords*>::iterator lit = ln.first; lit != ln.second; lit++) {
	 Funnel f(*lit, c, *rit);
	 std::cout << f.toString();
	 funnelQue.push_back(f);
	 rit++;
	 std::cout << "\n";
	 }
	std::cout << "\n starting while loop";

	while (!funnelQue.empty()) {
		std::cout << "fq.size = " << funnelQue.size() << std::endl;
		Funnel f = funnelQue.front();
		std::cout << "front" << std::endl;
		funnelQue.pop_front();
		std::cout << "popfront" << std::endl;
		std::pair<const Coords*, const Coords*> base = f.getBase();
		std::cout << "base" << std::endl;
		std::cout << base.first->toString() << "<->" << base.second->toString() << std::endl;
		const Coords* o = getOpposing(base.first, base.second, polygon);
		std::cout << "opposing" << o->toString() << std::endl;
		if (o != 0) {
			std::cout << "reacting" << std::endl;
			f.reactToOpposite(o, &funnelQue, neighbours);
			std::cout << "pushing back" << std::endl;
			funnelQue.push_back(f);
		} else {
			std::cout << "stopping" << std::endl;
		}
	}
}
std::list<const Coords*> findNeighbours(Coords* c) {
	std::list<const Coords*> neighbours { };
	std::pair<std::map<int, std::set<const Coords*>>::iterator, std::map<int, std::set<const Coords*>>::iterator> polyIt = c->getAllLeftN(); //Only keys area interesting here, first = begin iterator, second end iterator
	for (std::map<int, std::set<const Coords*>>::iterator it = polyIt.first; it != polyIt.second; it++) {
		findNeighboursInPolygon(c, it->first, &neighbours);
	}
	return neighbours;
}

int main() {

	//test();
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
				std::cout << "merginwith " << coordmap.find(c[i])->toString() << std::endl;
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

			//p.first->addToPolygon(0);
			cp[i]->addNeighbours(&*coordmap.find(c[l]), &*coordmap.find(c[r]), 0);
		}
	}
	std::cout << "coords found: " << coordmap.size() << std::endl;

	std::cout << "neighbours off 100 100" << std::endl;
	std::tr1::unordered_set<Coords>::iterator fit = coordmap.find(Coords(100, 100, 0));
	Coords* f = (Coords*) &*fit;
	std::pair<std::set<const Coords*>::iterator, std::set<const Coords*>::iterator> lneigh = f->getLeftNeighbours(0);
	std::pair<std::set<const Coords*>::iterator, std::set<const Coords*>::iterator> rneigh = f->getRightNeighbours(0);
	for (; lneigh.first != lneigh.second; lneigh.first++) {
		std::cout << "  " << (*lneigh.first)->toString() << "<===>" << (*rneigh.first)->toString() << std::endl;
		rneigh.first++;
	}
	std::tr1::unordered_set<Coords>::iterator lfit = coordmap.find(Coords(50, 200, 0));
	Coords* lf = (Coords*) &*lfit;
	std::tr1::unordered_set<Coords>::iterator rfit = coordmap.find(Coords(115, 160, 0));
	Coords* rf = (Coords*) &*rfit;

	std::cout << "opposing: " << getOpposing(lf, rf, 0)->toString();
	std::cout << "\ntest" << std::endl;
	std::cout << (Coords(300, 250) == Coords(300, 200)) << std::endl;

	std::cout << "neighbours:" << std::endl;
	std::list<const Coords*> n = findNeighbours(f);
	std::cout << "RESULT:" << std::endl;
	for (const Coords* c : n) {
		std::cout << c->toString() << std::endl;
	}

	return 0;
}

