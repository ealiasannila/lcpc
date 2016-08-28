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

double eucDistance(p2t::Point p1, p2t::Point p2) {
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

int addIntermidiatePoints(std::vector<p2t::Point*>* vec, std::vector<p2t::Point*>::iterator pit, std::vector<p2t::Point*>::iterator nextit, double maxDist) {
	p2t::Point p = **pit;
	p2t::Point next = **nextit;
	int pointsToAdd = std::ceil(eucDistance(p, next) / maxDist) - 1;
	std::vector<p2t::Point*> res(pointsToAdd);

	for (int i = 0; i < pointsToAdd; i++) {
		double x { ((next.x - p.x) / (pointsToAdd + 1)) * (i + 1) + p.x };
		double y { ((next.y - p.y) / (pointsToAdd + 1)) * (i + 1) + p.y };
		res[i] = new p2t::Point(x, y);
	}
	vec->insert(pit + 1, res.begin(), res.end());
	return pointsToAdd;
}

void intermidiatePoints(std::vector<std::vector<p2t::Point*>>* points, double maxDist) {
	for (unsigned int ring = 0; ring < points->size(); ring++) {
		for (unsigned int point = 0; point < points->at(ring).size(); point++) {
			// add intermidiate points if next point is too far.
			std::vector<p2t::Point*>::iterator next = (point + 1 != points->at(ring).size()) ? points->at(ring).begin() + point + 1 : points->at(ring).begin();
			if (eucDistance(*points->at(ring).at(point), **next) + 0.0000001 > maxDist) {
				addIntermidiatePoints(&points->at(ring), points->at(ring).begin() + point, next, maxDist);
			}
		}
	}
}

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
	finder.addPolygon(0,stp, p1, 1);
	finder.addPolygon(1,std::vector<p2t::Point*>{}, p2, 1);
	for (std::vector<p2t::Point*> vec : p1) {
		for (p2t::Point* p : vec) {
			delete p;
		}
	}

	std::cout << "HALFWAY\n";
	//finder.addPolygon(1, p2, 1, 0.5);

	std::vector<Coords> targets{Coords(0.5,2)};
	std::vector<Coords> path = finder.leastCostPath(Coords(0, 0), targets);
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

