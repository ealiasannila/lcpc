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
#include <iostream>
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

nContainer inline intersection(nContainer* ln, std::set<const Coords*>* rn) {
    nContainer res;
    for (auto li = ln->begin(); li != ln->end(); li++) {
        for (auto ri = rn->begin(); ri != rn->end(); ri++) {
            if (*ri>*li) {
                break;
            }
            if (*ri == *li) {
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
        double x{((next->x - p->x) / (pointsToAdd + 1)) * (i + 1) + p->x};
        double y{((next->y - p->y) / (pointsToAdd + 1)) * (i + 1) + p->y};
        res[i] = new p2t::Point(x, y);
    }
    vec->insert(pit + 1, res.begin(), res.end());
    return pointsToAdd;
}

void inline intermidiatePoints(std::vector<std::vector<p2t::Point*>>*points, double maxDist) {
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

int inline Orient(p2t::Point& pa, p2t::Point& pb, p2t::Point& pc) {
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

                point->x -= dx*ori;
                point->y -= dy*ori;
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
