//============================================================================
// Name        : lcpc.cpp
// Author      : Elias Annila
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "tests/test.h"
#include "../lib/geometry.h"
#include <algorithm>

Coords* getOpposing(Coords* l, Coords* r, int polygon) {
	std::cout<<"startOpposing";
	std::pair<std::set<Coords*>::iterator, std::set<Coords*>::iterator> ln =
			l->getLeftNeighbours(polygon);
	std::pair<std::set<Coords*>::iterator, std::set<Coords*>::iterator> rn =
			r->getRightNeighbours(polygon);
	for (std::set<Coords*>::iterator it = ln.first; it != ln.second; it++) {
		Coords* c = *it;
		std::cout << c->toString() << std::endl;
	}
	std::vector<Coords*> intersection;
	std::cout<<"startIntersection";
	std::set_intersection(ln.first, ln.second, rn.first, rn.second,
			std::back_inserter(intersection));
	std::cout<<"done" << std::endl;
	if (intersection.size() > 0) {
		return *intersection.begin();
	}

	std::cout<<"No opposing found"<<std::endl;
	return 0;
}

void findNeighboursInPolygon(Coords* c, int polygon,
		std::vector<Coords*>* neighbours) {
	std::pair<std::set<Coords*>::iterator, std::set<Coords*>::iterator> ln =
			c->getLeftNeighbours(polygon);
	std::pair<std::set<Coords*>::iterator, std::set<Coords*>::iterator> rn =
			c->getRightNeighbours(polygon);

	std::deque<Funnel> funnelQue { };
	std::set<Coords*>::iterator rit = rn.first;
	for (std::set<Coords*>::iterator lit = ln.first; lit != ln.second; lit++) {
		Funnel f(*lit, c, *rit);
		funnelQue.push_back(f);
		rit++;
	}
	while (!funnelQue.empty()) {
		std::cout << "fq.size = " << funnelQue.size() << std::endl;
		Funnel f = funnelQue.front();
		std::cout << "front" << std::endl;
		funnelQue.pop_front();
		std::cout << "popfront" << std::endl;
		std::pair<Coords*, Coords*> base = f.getBase();
		std::cout << "base" << std::endl;
		Coords* o = getOpposing(base.first, base.second, polygon);
		std::cout << "opposing" << std::endl;
		if (o != 0) {
			std::cout << "reacting" << std::endl;
			f.reactToOpposite(o, &funnelQue, neighbours);
			std::cout << "pushing back" << std::endl;
			funnelQue.push_back(f);
		}
	}
}
std::vector<Coords*> findNeighbours(Coords* c) {
	std::vector<Coords*> neighbours { };
	std::pair<std::map<int, std::set<Coords*>>::iterator,
			std::map<int, std::set<Coords*>>::iterator> polyIt =
			c->belongsToPolygons(); //Only keys area interesting here, first = begin iterator, second end iterator
	for (std::map<int, std::set<Coords*>>::iterator it = polyIt.first;
			it != polyIt.second; it++) {
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
	Triangles triangles = poly.triangles();

	Coords* start;

	for (std::list<Triangle>::iterator it = triangles.begin();
			it != triangles.end(); it++) {
		Triangle triangle = *it;
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
		}
		c[0].addNeighbours(&c[2], &c[1], 0);
		c[1].addNeighbours(&c[0], &c[2], 0);
		c[2].addNeighbours(&c[1], &c[0], 0);

		start = &c[1];
		std::cout<<start->toString<<std::endl;;
	}
	std::vector<Coords*> n = findNeighbours(start);
	for (Coords* c : n) {
		std::cout << c->toString() << std::endl;
	}

	return 0;
}

