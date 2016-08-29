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
	 *TODO Integrate to QGIS plugin
	 *TODO Implement A*
	 *TODO Writing results to shapefile (in QGIS)
	 *TODO Weakly simple polygon splitting
	 *TODO Performance...
	 *
	 */

	LcpFinder finder;

	std::vector<std::vector<p2t::Point*>> p1 = { { new p2t::Point { 0, 0 }, new p2t::Point(1, 0), new p2t::Point(1, 1), new p2t::Point(0, 1) } };
	std::vector<std::vector<p2t::Point*>> p2 = { { new p2t::Point { 0, 1 }, new p2t::Point(1, 1), new p2t::Point(1, 2), new p2t::Point(0, 2) } };
	std::cout << "vector created\n";
	std::cout << p1[0].size() << std::endl;
	intermidiatePoints(&p1, 0.5);
	intermidiatePoints(&p2, 0.5);
	std::cout << p1[0].size() << std::endl;
	std::cout << "intermidiate added\n";

	std::vector<p2t::Point*> stp = {new p2t::Point{0.5,0.5}};
	finder.addPolygon(p1, 1);
	std::cout<<"first added\n";
	finder.addPolygon(p2, 1);
	std::cout<<"second added\n";

	finder.addSteinerPoints(stp, 0);
	std::cout<<"Steiner added\n";
	std::tr1::unordered_set<Coords, CoordsHasher, std::equal_to<Coords>, std::allocator<Coords> > cmap = finder.getCoordmap();
	for(Coords c : cmap){
		std::cout<<c.toString()<<std::endl;
	}

	std::cout << "HALFWAY\n";
	//finder.addPolygon(1, p2, 1, 0.5);

	std::vector<Coords> targets{Coords(1,2)};
	std::vector<Coords> path = finder.leastCostPath(Coords(0.5, 0.5), targets);
	std::cout << "lcp done\n";
	for (Coords goal : path) {
		while (goal.getPred() != 0) {
			std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;
			goal = *goal.getPred();
		}
		std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;

	}

	return 0;
}

