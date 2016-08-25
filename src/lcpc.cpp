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
#include "lcpfinder.h"

long int line_id = 0;

class InsPoly: public Polygon {
private:
	Pointbase* prev = 0;
	Pointbase* first = 0;
public:
	void addPoint(unsigned int i, double x, double y);
	void addEdge(Pointbase* s, Pointbase* e);
	void closeRing();
};

void InsPoly::addPoint(unsigned int i, double x, double y) {
	Pointbase* point = new Pointbase(i, x, y, INPUT);
	if (x > _xmax)
		_xmax = x;
	if (x < _xmin)
		_xmin = x;
	if (y > _ymax)
		_ymax = y;
	if (y < _ymin)
		_ymin = y;
	_points[1] = point;

	if (first == 0) {
		first = point;
		prev = point;
	} else if (prev != 0) {
		Linebase* line = new Linebase(prev, point, INPUT);
		line_id++;
		_edges[line_id] = line;
		prev = point;
	}

}

void InsPoly::closeRing() {
	Linebase* line = new Linebase(prev, first, INPUT);
	line_id++;
	_edges[line_id] = line;

	prev = 0;
	first = 0;

}

int main() {
	/*
	 * HOW TO MAKE FAST:
	 * A*
	 * Data structures: Coordsin naapurimappi, Coordmappi arrayksi? Coords luokka kokonaan pois?
	 * Minheap, parempi update operaatio
	 *
	 */

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
	std::tr1::unordered_set<Coords>::iterator sit = coordmap.find(Coords(0, 0, 0));
	const Coords* s = &*sit;
	std::tr1::unordered_set<Coords>::iterator fit = coordmap.find(Coords(300, 250, 0));
	const Coords* f = &*fit;

	std::cout << "COORDMAP:" << std::endl;
	for (Coords c : coordmap) {
		std::cout << "cost: " << c.getToStart() << " pred: " << c.getPred() << std::endl;
	}
	return 0;
}

