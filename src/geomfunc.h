/*
 * geomfunc.h
 *
 *  Created on: Aug 28, 2016
 *      Author: elias
 */
#include "../lib/poly2tri.h"

#ifndef SRC_GEOMFUNC_H_
#define SRC_GEOMFUNC_H_

double inline eucDistance(const Coords* p1, const Coords* p2) {
	return sqrt(pow(p1->getX() - p2->getX(), 2) + pow(p1->getY() - p2->getY(), 2));
}

double inline eucDistance(p2t::Point p1, p2t::Point p2) {
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}



int inline addIntermidiatePoints(std::vector<p2t::Point*>* vec, std::vector<p2t::Point*>::iterator pit, std::vector<p2t::Point*>::iterator nextit, double maxDist) {
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

void  inline intermidiatePoints(std::vector<std::vector<p2t::Point*>>* points, double maxDist) {
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

#endif /* SRC_GEOMFUNC_H_ */
