/*
 * geomfunc.h
 *
 *  Created on: Aug 28, 2016
 *      Author: elias
 */

#ifndef SRC_GEOMFUNC_H_
#define SRC_GEOMFUNC_H_
#include "../lib/poly2tri.h"
#include "../lib/clipper/cpp/clipper.hpp"

#include <gdal/ogrsf_frmts.h>
#include "coords.h"
#include <iostream>
#include "defs.h"
#include "lcpfinder.h"

double inline eucDistance(const Coords* p1, const Coords* p2) {
    return std::sqrt(std::pow(p1->getX() - p2->getX(), 2) + std::pow(p1->getY() - p2->getY(), 2));
}

double inline eucDistance(p2t::Point* p1, p2t::Point* p2) {
    return std::sqrt(std::pow(p1->x - p2->x, 2) + std::pow(p1->y - p2->y, 2));
}

double inline eucDistance(p2t::Point* p1, const Coords* p2) {
    return std::sqrt(std::pow(p1->x - p2->getX(), 2) + std::pow(p1->y - p2->getY(), 2));
}


static const float A = 0.983398;
static const float B = 0.430664;
static const float C = 0.041;
static const float D = 0.08089;

float inline approxDistance(p2t::Point* p1, const Coords* p2) {

    float dx = p1->x - p2->getX();
    float dy = p1->y - p2->getY();
    float min, max, approx;
    if (dx < 0) dx = -dx;
    if (dy < 0) dy = -dy;

    if (dx < dy) {
        min = dx;
        max = dy;
    } else {
        min = dy;
        max = dx;
    }
    if (max > 4 * min) {
        return max;
    }
    return A * max + B * min - D * max;

}

nContainer inline intersection(std::vector<std::pair<const Coords*, double>>*ln, std::vector<std::pair<const Coords*, double>>*rn) {
    nContainer res;
    for (auto li = ln->begin(); li != ln->end(); li++) {
        for (auto ri = rn->begin(); ri != rn->end(); ri++) {
            /*
            if (*(ri->first)>*(li->first)) {
                break;
            }
             * */
            if (ri->first == li->first) {
                res.push_back(ri->first);
                if (res.size() == 2) {
                    return res;
                }
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
        double x{((next->x - p->x) / (pointsToAdd + 1)) * (i + 1) + p->x};
        double y{((next->y - p->y) / (pointsToAdd + 1)) * (i + 1) + p->y};
        res[i] = new p2t::Point(x, y);
    }
    vec->insert(pit + 1, res.begin(), res.end());

    return pointsToAdd;
}

int inline thirdCorner(const Triangle* t, const Coords* c1, const Coords* c2) {
    for (int i = 0; i < 3; i++) {
        if (t->points[i] != c1 and t->points[i] != c2) {

            return i;
        }
    }
    return -1;
}

std::pair<int, const Triangle *> inline commonTriangle(const Coords* c1, const Coords* c2, const Triangle* currentTriangle) {
    for (int nextPolygon : c1->belongsToPolygons()) {
        if (c2->belongsToPolygon(nextPolygon)) {
            for (auto it = c1->getTriangles(nextPolygon)->begin(); it != c1->getTriangles(nextPolygon)->end(); it++) {
                const Triangle * t1 = *it;
                if (t1 == currentTriangle) {
                    continue;
                }
                auto it2 = std::find(c2->getTriangles(nextPolygon)->begin(), c2->getTriangles(nextPolygon)->end(), t1);
                if (it2 != c2->getTriangles(nextPolygon)->end()) {

                    return std::make_pair(nextPolygon, t1);
                }
            }
        }
    }
}

void inline intermidiatePointsForRing(std::vector<p2t::Point*>*points, double maxDist, bool isring) {
    for (unsigned int point = 0; point < points->size(); point++) {
        // add intermidiate points if next point is too far.
        if (!isring and point + 1 == points->size()) {
            break;
        }
        std::vector<p2t::Point*>::iterator next = (point + 1 != points->size()) ? points->begin() + point + 1 : points->begin();
        if (eucDistance(points->at(point), *next) + 0.0000001 > maxDist) {

            addIntermidiatePoints(points, points->begin() + point, next, maxDist);
        }
    }
}

void inline intermidiatePoints(std::vector<std::vector < p2t::Point*>>*points, double maxDist, bool isring) {
    for (unsigned int ring = 0; ring < points->size(); ring++) {

        intermidiatePointsForRing(&points->at(ring), maxDist, isring);
    }
}

void inline intermidiatePoints(std::vector<std::vector < p2t::Point*>>*points, double maxDist) {

    return intermidiatePoints(points, maxDist, true);
}

bool inline compDijkstra(const Coords* x, const Coords * y) {

    return (x->getToStart() > y->getToStart());
}

bool inline compAstar(const Coords* x, const Coords * y) {

    return ((x->getToStart() + x->getToEnd()) > (y->getToStart() + y->getToEnd()));

}

bool inline coordsInTriangle(Triangle* t, const Coords * c) {
    for (int i = 0; i < 3; i++) {
        if (c->isRight(t->points[i], t->points[(i + 1) % 3]) == 1) {

            return false;
        }
    }
    return true;
}

/*
 * 0 = no crossing
 * 1 = crossing
 * 2 = crossing at endpoint
 */
int inline segmentCrossing(const Coords* a1, const Coords* a2, const Coords* b1, const Coords* b2) {
    int a1Tob = a1->isRight(b1, b2);
    int a2Tob = a2->isRight(b1, b2);
    int b1Toa = b1->isRight(a1, a2);
    int b2Toa = b2->isRight(a1, a2);
    if (a1Tob == a2Tob or b1Toa == b2Toa) {
        return 0;
    }
    if (a1Tob == -1 * a2Tob and b1Toa == -1 * b2Toa) {
        return 1;
    }
    return 2;
}


void inline insertToNset(nSet* nset, const Coords* a, double fric, const Coords* b) {
    double cost {b->getToStart() + eucDistance(b, a) * fric};
    auto p = nset->insert(std::make_pair(a,cost));
    if (!p.second and p.first->second > cost) {
        p.first->second = cost;
    }
}


bool inline inside(std::vector<std::vector < p2t::Point*>> polygon, p2t::Point * point) {
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

int inline Orient(p2t::Point& pa, p2t::Point& pb, p2t::Point & pc) {
    double detleft = (pa.x - pc.x) * (pb.y - pc.y);
    double detright = (pa.y - pc.y) * (pb.x - pc.x);
    double val = detleft - detright;
    if (val > -0.000001 && val < 0.000001) {
        return 0;
    } else if (val > 0) {

        return -1;
    }
    return 1;
}

bool inline pointOnSegment(p2t::Point* l1, p2t::Point* l2, p2t::Point * point) {
    double cross = (point->y - l1->y) * (l2->x - l1->x) - (point->x - l1->x) * (l2->y - l1->y);
    if (std::abs(cross) > 0.000001) {
        return false;
    }
    double dot = (point->x - l1->x) * (l2->x - l1->x) + (point->y - l1->y)*(l2->y - l1->y);
    if (dot < 0) {
        return false;
    }
    double sqrd = (l2->x - l1->x)*(l2->x - l1->x) + (l2->y - l1->y)*(l2->y - l1->y);
    if (dot > sqrd) {

        return false;
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
    std::map<std::pair<int, int>, int> points;

    std::vector<unsigned int> holes;
    std::vector<unsigned int> outers;
    for (int i = 0; i < paths.size(); i++) {
        ClipperLib::Path path = paths[i];
        if (ClipperLib::Orientation(path)) {
            std::vector < p2t::Point*> outer;
            for (ClipperLib::IntPoint ip : path) {
                if (!points.insert(std::make_pair(std::make_pair(ip.X, ip.Y), 0)).second) {
                    std::cout << "DOUBLE POINT ON OUTER\n";
                }
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
            auto ins = points.insert(std::make_pair(std::make_pair(ip.X, ip.Y), hole));

            p2t::Point* point = new p2t::Point{(double) ip.X / scale, (double) ip.Y / scale};
            if (!ins.second) {
                std::cout << "DOUBLE POINT ON RING " << ins.first->second << "\n";
                std::cout << "Current ring: " << hole << std::endl;
                std::cout << "xy: " << point->x << "," << point->y << std::endl;
                int j = 1;
                int ori = 0;
                p2t::Point next;
                p2t::Point prev;
                do {
                    int pipInd = paths[hole].size() - 1 - i + j;

                    int nipInd = paths[hole].size() - 1 - i - j;
                    if (nipInd < 0) {
                        nipInd = paths[hole].size() - nipInd;
                    }
                    if (pipInd > paths[hole].size() - 1) {
                        pipInd = pipInd - paths[hole].size();
                    }
                    j++;
                    ClipperLib::IntPoint nip = paths[hole].at(nipInd);
                    ClipperLib::IntPoint pip = paths[hole].at(pipInd);
                    next = p2t::Point{(double) nip.X / scale, (double) nip.Y / scale};
                    prev = p2t::Point{(double) pip.X / scale, (double) pip.Y / scale};
                    ori = Orient(prev, *point, next);
                } while (ori == 0);


                p2t::Point moved{(next.x - prev.x) / 2 + prev.x, (next.y - prev.y) / 2 + prev.y};
                double e = eucDistance(point, &moved);
                double dx = (moved.x - point->x) / e;
                double dy = (moved.y - point->y) / e;

                point->x -= dx * ori;
                point->y -= dy * ori;
            }
            inner.push_back(point);
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



std::array<double, 2> inline lineIntersection(const Coords* l1, const Coords* l2, const Coords* b1, const Coords * b2) {
    double slopeL = (l2->getY() - l1->getY()) / (l2->getX() - l1->getX());
    double slopeB = (b2->getY() - b1->getY()) / (b2->getX() - b1->getX());
    double interceptL = l1->getY() - slopeL * l1->getX();
    double interceptB = b1->getY() - slopeB * b1->getX();
    double x = (interceptL - interceptB) / (slopeB - slopeL);
    double y = slopeB * x + interceptB;

    if (l1->getX() == l2->getX()) {
        x = l1->getX();
        if (b1->getY() != b2->getY()) {
            y = slopeB * x + interceptB;
        }
    } else if (l1->getY() == l2->getY()) {
        y = l1->getY();
        if (b1->getX() != b2->getX()) {
            x = (y - interceptB) / slopeB;
        }
    }
    if (b1->getX() == b2->getX()) {
        x = l1->getX();
        if (l1->getY() != l2->getY()) {
            y = slopeL * x + interceptL;
        }
    } else if (b1->getY() == b2->getY()) {
        y = b1->getY();
        if (l1->getX() != l2->getX()) {

            x = (y - interceptL) / slopeL;
        }
    }
    return std::array<double, 2>{x, y};
}

std::vector<std::vector<std::vector < p2t::Point*>>> inline dumbSimplify(OGRPolygon * polygon) {
    std::vector<std::vector<std::vector < p2t::Point*>>> out;
    out.push_back(std::vector<std::vector < p2t::Point*>>
    {
    });
    OGRLineString* ext = polygon->getExteriorRing();
    out.back().push_back(std::vector<p2t::Point*>{});

    for (int i = 0; i < ext->getNumPoints(); i++) {
        int index = ext->getNumPoints() - 1 - i;
        out.back().back().push_back(new p2t::Point{ext->getX(i), ext->getY(i)});
    }
    p2t::Point* d = out.back().back().back();
    p2t::Point* s = out.back().back()[0];
    if (s->x == d->x and s->y == d->y) {
        out.back().back().pop_back();
        delete d;
    }

    for (int i = 0; i < polygon->getNumInteriorRings(); i++) {
        out.back().push_back(std::vector<p2t::Point*>{});
        OGRLineString* inter = polygon->getInteriorRing(i);
        for (int j = 0; j < inter->getNumPoints(); j++) {
            out.back().back().push_back(new p2t::Point{inter->getX(j), inter->getY(j)});
        }
        p2t::Point* d = out.back().back().back();
        p2t::Point* s = out.back().back()[0];
        if (s->x == d->x and s->y == d->y) {
            out.back().back().pop_back();
            delete d;
        }
    }
    return out;

}



#endif /* SRC_GEOMFUNC_H_ */
