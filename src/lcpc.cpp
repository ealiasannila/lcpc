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

int main() {
	/*
	 * HOW TO MAKE FAST:
	 * A*
	 * Data structures: Coordsin naapurimappi, Coordmappi arrayksi? Coords luokka kokonaan pois?
	 * Minheap, parempi update operaatio
	 *
	 */

	std::cout << "Starting\n";
	Polygon poly { };
	std::cout << "Constructed\n";
	poly.setDebugOption(false);      //set debug flag;
	std::cout << "Debug = false\n";
	poly.addPoint(1, 0, 0, 0);
	std::cout << "P1\n";
	poly.addPoint(2, 1, 0, 0);
	std::cout << "P2\n";
	poly.addPoint(3, 1, 1, 0);
	std::cout << "P3\n";
	poly.addPoint(4, 1.5, 0.5, 0);
	std::cout << "P4\n";
	poly.addPoint(5, 0, 1, 0);
	std::cout << "P5\n";
	poly.addEdges();
	std::cout << "EDGES\n";

	for(std::pair<int, Linebase*> pair : poly.edges()){
		std::cout<<"key: "<<pair.first<<std::endl;
		std::cout<<pair.second->endPoint(0)->id<<"--"<<pair.second->endPoint(1)->id<<std::endl;
	}


	std::cout<<"prevtest:"<<std::endl;
	std::cout<<poly.prev(1)<<std::endl;
	std::cout<<poly.prev(2)<<std::endl;
	std::cout<<poly.prev(3)<<std::endl;
	std::cout<<poly.prev(4)<<std::endl;
	std::cout<<poly.prev(5)<<std::endl;

	std::cout<<"nexttest:"<<std::endl;
	std::cout<<poly.next(1)<<std::endl;
	std::cout<<poly.next(2)<<std::endl;
	std::cout<<poly.next(3)<<std::endl;
	std::cout<<poly.next(4)<<std::endl;
	std::cout<<poly.next(5)<<std::endl;

	poly.initializate();
	std::cout << "INIT\n";

	poly.triangulation();
	std::cout << "triangulaton done" << std::endl;
	PointbaseMap points = poly.points();
	std::tr1::unordered_set<Coords, CoordsHasher> coordmap;

	Triangles triangles = poly.triangles();

	int t1 = 1;
	for (std::list<Triangle>::iterator it = triangles.begin(); it != triangles.end(); it++) {
		std::cout << "\n\n--Triangle: " << t1 << "--\n";
		Triangle triangle = *it;
		for(unsigned int i : triangle){
			std::cout<<i<<std::endl;
		}
	}
	return 0;

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
	/*
	 //set start vertex:
	 std::tr1::unordered_set<Coords>::iterator sit = coordmap.find(Coords(0, 0, 0));
	 const Coords* s = &*sit;
	 std::tr1::unordered_set<Coords>::iterator fit = coordmap.find(Coords(300, 250, 0));
	 const Coords* f = &*fit;

	 std::cout << "COORDMAP:" << std::endl;
	 for (Coords c : coordmap) {
	 std::cout << "cost: " << c.getToStart() << " pred: " << c.getPred() << std::endl;
	 }
	 */
	return 0;
}

