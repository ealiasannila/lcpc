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
	 *
	 *TODO Add maximum node distance
	 *TODO Add points inside polygons (as stainer points?)
	 *TODO Integrate to QGIS plugin
	 *TODO Implement A*
	 *TODO Writing results to shapefile (in QGIS)
	 *TODO Weakly simple polygon splitting
	 *TODO Performance...
	 *
	 */

	LcpFinder finder;

	std::vector<std::vector<Coords>> p1 = { { Coords(0, 0), Coords(1, 0), Coords(1, 1), Coords(0, 1) } };
	std::vector<std::vector<Coords>> p2 = { { Coords(0, 1), Coords(1, 1), Coords(1, 2), Coords(0, 2) } };
	finder.addPolygon(0, p1, 1, 0.5);
	finder.addPolygon(1, p2, 1, 0.5);
	std::vector<Coords> path = finder.leastCostPath(Coords(0, 0), Coords(1, 2));
	std::cout<<"lcp done\n";
	for (Coords goal : path) {
		while (goal.getPred() != 0) {
			std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;
			goal = *goal.getPred();
		}
		std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;

	}

	return 0;
}

