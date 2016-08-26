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
#include "lcpfinder.h"

int main() {
	/*
	 * HOW TO MAKE FAST:
	 * A*
	 * Data structures: Coordsin naapurimappi, Coordmappi arrayksi? Coords luokka kokonaan pois?
	 * Minheap, parempi update operaatio
	 *
	 */

	LcpFinder finder;

	std::vector<std::vector<Coords>> p1 = { { Coords(0, 0), Coords(1, 0), Coords(1, 1), Coords(0, 1) } };
	std::vector<std::vector<Coords>> p2 = { { Coords(1, 0), Coords(2, 0), Coords(2, 1), Coords(1, 1) } };
	finder.addPolygon(0, p1);
	finder.addPolygon(1, p2);
		std::vector<Coords> path = finder.leastCostPath(Coords(0, 0), Coords(2, 1));
	for (Coords goal : path) {
		while (goal.getPred() != 0) {
			std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;
			goal = *goal.getPred();
		}
		std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;

	}

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

