/*
 * geomfunc.h
 *
 *  Created on: Aug 28, 2016
 *      Author: elias
 */
#include "../lib/poly2tri.h"
#include "../lib/clipper/cpp/clipper.hpp"

#include <gdal/ogrsf_frmts.h>
#include "coords.h"

#ifndef SRC_GEOMFUNC_H_
#define SRC_GEOMFUNC_H_

double inline eucDistance(const Coords* p1, const Coords* p2) {
	return std::sqrt(std::pow(p1->getX() - p2->getX(), 2) + std::pow(p1->getY() - p2->getY(), 2));
}

double inline eucDistance(p2t::Point* p1, p2t::Point* p2) {
	return std::sqrt(std::pow(p1->x - p2->x, 2) + std::pow(p1->y - p2->y, 2));
}


double inline eucDistance(p2t::Point* p1, const Coords* p2) {
	return std::sqrt(std::pow(p1->x - p2->getX(), 2) + std::pow(p1->y - p2->getY(), 2));
}


nContainer inline intersection(nContainer* ln, SortedVector<const Coords*>* rn){
    nContainer res;
    for(auto li = ln->begin(); li != ln->end(); li++){
        for (auto ri = rn->begin(); ri != rn->end(); ri++){
            if(*ri>*li){
                break;
            }
            if(*ri==*li){
                res.push_back(*ri);
            }
        }
    }
    return res;
}

int inline addIntermidiatePoints(std::vector<p2t::Point*>* vec, std::vector<p2t::Point*>::iterator pit, std::vector<p2t::Point*>::iterator nextit, double maxDist) {
	p2t::Point* p = *pit;
	p2t::Point* next = *nextit;
	int pointsToAdd = std::ceil(eucDistance(p, next) / maxDist) - 1;
	std::vector<p2t::Point*> res(pointsToAdd);

	for (int i = 0; i < pointsToAdd; i++) {
		double x { ((next->x - p->x) / (pointsToAdd + 1)) * (i + 1) + p->x };
		double y { ((next->y - p->y) / (pointsToAdd + 1)) * (i + 1) + p->y };
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
			if (eucDistance(points->at(ring).at(point), *next) + 0.0000001 > maxDist) {
				addIntermidiatePoints(&points->at(ring), points->at(ring).begin() + point, next, maxDist);
			}
		}
	}
}

bool inline compDijkstra(const Coords* x, const Coords* y) {
    return (x->getToStart() > y->getToStart());
}

bool inline compAstar(const Coords* x, const Coords* y) {
    return ((x->getToStart() + x->getToEnd()) > (y->getToStart() + y->getToEnd()));

}

bool inline inside(std::vector<std::vector<p2t::Point*>> polygon, p2t::Point* point) {
    bool exterior = true;
    for (std::vector<p2t::Point*> ring : polygon) {

        int i, j = 0;
        bool c = false;
        for (i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
            p2t::Point* ringPoint = ring[i];
            p2t::Point* ringPoint2 = ring[j];
            if (((ringPoint->y > point->y) != (ringPoint2->y > point->y)) &&
                    (point->x < (ringPoint2->x - ringPoint->x) * (point->y - ringPoint->y) / (ringPoint2->y - ringPoint->y) + ringPoint->x))
                c = !c;
        }
        if (exterior and !c) {
            return false;
        }
        if (!exterior and c) {
            return false;
        }
        exterior = false;
    }
    return true;
}

std::vector<std::vector<std::vector < p2t::Point*>>> inline simplify(OGRPolygon * polygon) {
    int scale = 10000;
    ClipperLib::Paths paths;
    paths.push_back(ClipperLib::Path{});
    OGRLineString* ext = polygon->getExteriorRing();
    for (int i = 0; i < ext->getNumPoints(); i++) {
        int index = ext->getNumPoints() - 1 - i;
        paths.back().push_back(ClipperLib::IntPoint(ext->getX(index) * scale, ext->getY(index) * scale));
    }


    for (int i = 0; i < polygon->getNumInteriorRings(); i++) {
        paths.push_back(ClipperLib::Path{});
        OGRLineString* inter = polygon->getInteriorRing(i);
        for (int j = 0; j < inter->getNumPoints(); j++) {
            paths.back().push_back(ClipperLib::IntPoint(inter->getX(j) * scale, inter->getY(j) * scale));
        }
    }

    std::vector<std::vector<std::vector < p2t::Point*>>> out;


    ClipperLib::SimplifyPolygons(paths, ClipperLib::pftEvenOdd);

    std::vector<unsigned int> holes;
    std::vector<unsigned int> outers;
    for (int i = 0; i < paths.size(); i++) {
        ClipperLib::Path path = paths[i];
        if (ClipperLib::Orientation(path)) {
            std::vector < p2t::Point*> outer;
            for (ClipperLib::IntPoint ip : path) {
                outer.push_back(new p2t::Point{(double) ip.X / scale, (double) ip.Y / scale});
            }
            out.push_back(std::vector<std::vector < p2t::Point*>>
            {
                outer
            });
        } else {
            holes.push_back(i);
        }
    }

    for (unsigned int hole : holes) {
        std::vector < p2t::Point*> inner;

        for (unsigned int i = 0; i < paths[hole].size(); i++) {
            ClipperLib::IntPoint ip = paths[hole].at(paths[hole].size() - 1 - i);
            inner.push_back(new p2t::Point{(double) ip.X / scale, (double) ip.Y / scale});
        }
        unsigned int outIndex = 0;
        for (std::vector<std::vector < p2t::Point*>> outer : out) {
            if (inside(outer, inner[0])) {
                out.at(outIndex).push_back(inner);
            }
            outIndex++;
        }
    }
    return out;

}


#endif /* SRC_GEOMFUNC_H_ */
