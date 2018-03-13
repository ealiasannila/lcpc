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

#define START 0
#define STOP 1
#define SPLIT 2
#define MERGE 3
#define ORDINARY 4

double inline eucDistance(const Coords* p1, const Coords* p2) {
    return std::sqrt(std::pow(p1->getX() - p2->getX(), 2) + std::pow(p1->getY() - p2->getY(), 2));
}

double inline eucDistance(std::array<double, 2> p1, std::array<double, 2> p2) {
    return std::sqrt(std::pow(p1[0] - p2[0], 2) + std::pow(p1[1] - p2[1], 2));
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

int inline neighbouringPolygon(const Coords* a, const Coords* b, int polygon) {
    std::vector<int> ap = a->belongsToPolygons();
    std::vector<int> bp = b->belongsToPolygons();
    for (int i : ap) {
        for (int j : bp) {
            if (j == i and j != polygon) {
                return j;
            }
        }
    }
    return -1;
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

void inline insertToNset(nSet* nset, const Coords* a, double fric, const Coords* b, double crossingCost) {
    double cost{b->getToStart() + eucDistance(b, a) * fric + crossingCost};
    auto p = nset->insert(std::make_pair(a, cost));
    if (!p.second and p.first->second > cost) {
        p.first->second = cost;
    }

}

void inline insertToNset(nSet* nset, const Coords* a, double fric, const Coords* b) {
    insertToNset(nset, a, fric, b, 0);
}

bool inline inside(std::vector<std::vector <const Coords*>> polygon, p2t::Point * point) {
    bool exterior = true;
    for (std::vector<const Coords*> ring : polygon) {

        int i, j = 0;
        bool c = false;
        for (i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
            const Coords* ringPoint = ring[i];
            const Coords* ringPoint2 = ring[j];
            if (((ringPoint->getY() > point->y) != (ringPoint2->getY() > point->y)) &&
                    (point->x < (ringPoint2->getX() - ringPoint->getX()) * (point->y - ringPoint->getY()) / (ringPoint2->getY() - ringPoint->getY()) + ringPoint->getX()))
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

bool inline inside(std::vector<std::vector < p2t::Point>> polygon, p2t::Point point) {
    bool exterior = true;
    for (std::vector<p2t::Point> ring : polygon) {

        int i, j = 0;
        bool c = false;
        for (i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
            p2t::Point ringPoint = ring[i];
            p2t::Point ringPoint2 = ring[j];
            if (((ringPoint.y > point.y) != (ringPoint2.y > point.y)) &&
                    (point.x < (ringPoint2.x - ringPoint.x) * (point.y - ringPoint.y) / (ringPoint2.y - ringPoint.y) + ringPoint.x))
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

std::vector<std::vector<std::vector < p2t::Point>>> inline simplify(OGRPolygon * polygon) {
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

    std::vector<std::vector<std::vector < p2t::Point>>> out;


    ClipperLib::SimplifyPolygons(paths, ClipperLib::pftEvenOdd);
    std::map<std::pair<int, int>, int> points;

    std::vector<unsigned int> holes;
    std::vector<unsigned int> outers;
    for (int i = 0; i < paths.size(); i++) {
        ClipperLib::Path path = paths[i];
        if (ClipperLib::Orientation(path)) {
            std::vector < p2t::Point> outer;
            for (ClipperLib::IntPoint ip : path) {
                if (!points.insert(std::make_pair(std::make_pair(ip.X, ip.Y), 0)).second) {
                    std::cout << "DOUBLE POINT ON OUTER\n";
                }
                outer.push_back(p2t::Point{(double) ip.X / scale, (double) ip.Y / scale});
            }
            out.push_back(std::vector<std::vector < p2t::Point>>
            {
                outer
            });
        } else {
            holes.push_back(i);
        }
    }

    for (unsigned int hole : holes) {
        std::vector < p2t::Point> inner;
        for (unsigned int i = 0; i < paths[hole].size(); i++) {
            ClipperLib::IntPoint ip = paths[hole].at(paths[hole].size() - 1 - i);
            auto ins = points.insert(std::make_pair(std::make_pair(ip.X, ip.Y), hole));

            p2t::Point point{(double) ip.X / scale, (double) ip.Y / scale};
            if (!ins.second) {
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
                    ori = Orient(prev, point, next);
                } while (ori == 0);


                p2t::Point moved{(next.x - prev.x) / 2 + prev.x, (next.y - prev.y) / 2 + prev.y};
                double e = eucDistance(&point, &moved);
                double dx = (moved.x - point.x) / e;
                double dy = (moved.y - point.y) / e;

                point.x -= dx * ori;
                point.y -= dy * ori;
            }
            inner.push_back(point);
        }
        unsigned int outIndex = 0;
        for (std::vector<std::vector < p2t::Point >> outer : out) {
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

void inline printTriangle2(const Triangle* t) {
    std::cout << std::fixed;
    std::cout << "Triangle : " << t << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "p: " << t->points[i]->getX() << "," << t->points[i]->getY() << std::endl;
        std::cout << "n: " << t->neighbours[i] << std::endl;
    }
}

bool inline cmpCoordsPointersByXY(const Coords* a, const Coords* b) {
    return (a->getX() < b->getX()) or(a->getX() == b->getX() and a->getY() < b->getY());
}

int inline getType(const Coords* c, std::unordered_map<const Coords*, std::vector<const Coords*>>*neighbors) {
    const Coords* l = neighbors->at(c)[0];
    const Coords* r = neighbors->at(c)[1];

    bool compR = cmpCoordsPointersByXY(r, c);
    bool compL = cmpCoordsPointersByXY(l, c);

    if (compL and compR) {
        if (r->isRight(l, c) == -1) {
            return MERGE;
        }
        return STOP;
    }
    if (!compL and !compR) {
        if (r->isRight(l, c) == -1) {
            return SPLIT;
        }
        return START;
    }
    return ORDINARY;


}

bool inline cmp(Edge a, Edge b) {
    if (a.first->isRight(b.first, b.second) == 0 and a.second->isRight(b.first, b.second) == 0) {
        return false;
    }

    int o1 = a.first->isRight(b.first, b.second);
    int o2 = a.second->isRight(b.first, b.second);
    if (o1 != -1 and o2 != -1) {
        return true;
    }
    int o3 = b.first->isRight(a.first, a.second);
    int o4 = b.second->isRight(a.first, a.second);
    return o3 != 1 and o4 != 1;

}

auto inline createEdge(const Coords* a, const Coords * b) {
    if (a == b or a == nullptr or b == nullptr) {
        std::cout << "FUCKED UP" << std::endl;
        exit(1);
    }
    if (cmpCoordsPointersByXY(a, b)) {
        return std::make_pair(a, b);
    }
    return std::make_pair(b, a);
}

auto inline setItDown(Edge lower, std::map<Edge, const Coords*, bool (*)(Edge, Edge)>*edges) {
    auto it = edges->lower_bound(lower);

    if (it != edges->end() and it != edges->begin()) {
        it--;
    } else if (it == edges->begin()) {
        it = edges->end();
    }
    return it;
}

auto inline leftMost(const Coords* c, const Coords * prev, std::unordered_map<const Coords*, std::vector<const Coords*>>*neighbors) {
    double maxAngle = std::numeric_limits<double>::min();
    const Coords* leftMost = nullptr;

    for (const Coords* n : neighbors->at(c)) {
        if (n == prev) {
            continue;
        }
        double x1 = c->getX() - prev->getX();
        double y1 = c->getY() - prev->getY();
        double x2 = c->getX() - n->getX();
        double y2 = c->getY() - n->getY();

        double dot = x1 * x2 + y1*y2;
        double det = x1 * y2 - y1*x2;
        double angle = atan2(det, dot);
        if (angle < 0) {
            angle = 2 * M_PI + angle;
        }
        if (angle > maxAngle) {
            maxAngle = angle;
            leftMost = n;
        }

    }
    if (leftMost == nullptr) {
        std::cout << "leftmost not found" << std::endl;
        exit(1);
    }
    return leftMost;
}

auto inline addDiagonal(Edge e, std::vector<Edge>*diagonals, std::vector<Edge>*starts, std::unordered_map<const Coords*, std::vector<const Coords*>>*neighbors) {
    int type = getType(e.first, neighbors);

    if (type != START and type != MERGE) {
        auto leftMostNeighbor = leftMost(e.first, e.second, neighbors);
        if (cmpCoordsPointersByXY(leftMostNeighbor, e.first)) {
            starts->push_back(e);
        } else {
            starts->push_back(createEdge(e.first, leftMostNeighbor));
        }

    }
    neighbors->at(e.first).push_back(e.second);
    neighbors->at(e.second).push_back(e.first);

    diagonals->push_back(e);


};

auto inline splitToMonotone(std::vector<std::vector<const Coords*>> polygon, std::vector<Edge>* diagonals) {

    std::vector<const Coords*> combined{};
    std::unordered_map<const Coords*, std::vector<const Coords*>> neighbors;


    for (auto ring = polygon.begin(); ring != polygon.end(); ring++) {
        const Coords* prev2 = ring->at(ring->size() - 2);
        const Coords* prev1 = ring->at(ring->size() - 1);
        for (auto c = ring->begin(); c != ring->end(); c++) {
            neighbors.insert(std::make_pair(prev1, std::vector<const Coords*>{*c, prev2}));
            combined.push_back(*c);
            prev2 = prev1;
            prev1 = *c;
        }
    }


    std::sort(combined.begin(), combined.end(), cmpCoordsPointersByXY);
    std::map<Edge, const Coords*, bool (*)(Edge, Edge) > edges(cmp);
    std::vector<Edge> starts{};




    
    for (const Coords* c : combined) {
        int vtype = getType(c, &neighbors);


        if (vtype == START) {
            starts.push_back(createEdge(c, neighbors.at(c)[0]));
        }

        //IF vtype in split, merge, start or stop: higher = higher, if vtype is ordinary higher = righmost edge
        Edge higher = createEdge(c, neighbors.at(c)[0]);
        Edge lower = createEdge(c, neighbors.at(c)[1]);
        bool orient = *(neighbors.at(c)[0]) < *(neighbors.at(c)[1]);

        auto itUp = edges.end();
        auto itDown = edges.end();






        /*
         Messed up section setting itUp and Down
         */
        if (vtype != ORDINARY) {
            if (cmp(higher, lower)) {
                Edge tmp = lower;
                lower = higher;
                higher = tmp;

            }
            if (vtype == SPLIT or vtype == MERGE) {
                itUp = edges.upper_bound(higher);
                itDown = setItDown(lower, &edges);
            }
        } else if (orient) {
            Edge tmp = lower;
            lower = higher;
            higher = tmp;
            itDown = setItDown(lower, &edges);
        } else {
            itUp = edges.upper_bound(lower);

        }
        

        //OK section setting diagonals
        if (vtype == SPLIT) {
            if (itUp == edges.end()) {
                std::cout << "mysterious split.." << std::endl;
                exit(1);
            } else {
                addDiagonal(createEdge(c, itUp->second), diagonals, &starts, &neighbors);
            }
        }
        if (itDown != edges.end()) {
            if (getType(itDown->second, &neighbors) == MERGE) {
                addDiagonal(createEdge(c, itDown->second), diagonals, &starts, &neighbors);
            }
            edges[itDown->first] = c;
        }

        if (itUp != edges.end()) {
            if (getType(itUp->second, &neighbors) == MERGE) {
                addDiagonal(createEdge(c, itUp->second), diagonals, &starts, &neighbors);
            }
            edges[itUp->first] = c;
        }




        //OK section updating edges
        if (vtype == SPLIT or vtype == START) {

            edges[higher] = c;
            edges[lower] = c;
        } else if (vtype == STOP or vtype == MERGE) {

            edges.erase(higher);
            edges.erase(lower);
        } else if (vtype == ORDINARY) {
            edges.erase(lower);
            edges[higher] = c;
        }

   
    }
   
    std::vector<std::vector<const Coords*>> monotone
    {
    };
    monotone.reserve(starts.size());

    for (Edge e : starts) {
        //matplotlibcpp::annotate("S", e.first->getX(), e.first->getY());
        const Coords* prev = e.first;
        const Coords* s = prev;
        const Coords* c = e.second;
        std::vector<const Coords*> mono{s};

        while (c != s) {
            mono.push_back(c);
            const Coords* next = leftMost(c, prev, &neighbors);
            prev = c;
            c = next;
        }

        monotone.push_back(mono);
    }
    return monotone;
}

auto inline triangulateMonotone(std::vector<const Coords*> monotone, int polygonId) {

    //Neighbor information could be retrieved from splitToMonotone operation at O(1)
    std::vector<std::pair<const Coords*, int>> sorted;
    std::map<const Coords*, const Coords*> nextAlongMonotone{};
    sorted.reserve(monotone.size());
    for (int i = 0; i < monotone.size(); i++) {
        auto c = monotone[i];
        auto r = monotone[(i - 1 + monotone.size()) % monotone.size()];
        auto l = monotone[(i + 1) % monotone.size()];
        nextAlongMonotone[c] = r;
        if (*l<*c and *c<*r) {
            sorted.push_back(std::make_pair(c, -1));

        } else if (*l>*c and *c>*r) {
            sorted.push_back(std::make_pair(c, 1));

        } else {
            sorted.push_back(std::make_pair(c, 0));

        }
    }

    auto cmpCoordsPair = [](std::pair<const Coords*, int>a, std::pair<const Coords*, int> b) {
        return cmpCoordsPointersByXY(a.first, b.first);
    };
    std::sort(sorted.begin(), sorted.end(), cmpCoordsPair);
    std::vector<std::pair<const Coords*, int>> stack
    {
        sorted[0], sorted[1]
    };


    std::vector<Triangle*> triangles{};

    std::set<Triangle*> triangleSet{};
    auto targetNumberOfNeighbors = [&nextAlongMonotone](Triangle * t) {
        int target = 3;
        for (int i = 0; i < 3; i++) {
            if (nextAlongMonotone[t->points[i]] == t->points[(i + 2) % 3]) {
                target--;
            }
        }
        return target;
    };

    auto trianglesAreNeighbors = [](Triangle* t1, Triangle * t2) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (t1->points[i] == t2->points[j] and t1->points[(i + 1) % 3] == t2->points[(j + 2) % 3]) {

                    t1->neighbours[(i + 2) % 3] = t2;
                    t2->neighbours[(j + 1) % 3] = t1;
                    return true;
                }
            }
        }
        return false;
    };


    auto neighborsFull = [targetNumberOfNeighbors](Triangle * t) {
        int missing = 0;
        for (auto n : t->neighbours) {
            if (n == nullptr) {
                missing++;
            }
        }
        return targetNumberOfNeighbors(t) + missing == 3;
    };

    auto addToTriangles = [&targetNumberOfNeighbors,&polygonId, &triangles, &triangleSet, &trianglesAreNeighbors, &neighborsFull](std::pair<const Coords*, int>c, std::pair<const Coords*, int>k, std::pair<const Coords*, int>st) {
        triangles.push_back(new Triangle(std::array<const Coords*, 3>{k.first, c.first, st.first}));
        c.first->addTriangle(triangles.back(), polygonId);
        k.first->addTriangle(triangles.back(), polygonId);
        st.first->addTriangle(triangles.back(), polygonId);
        int neighbors = 0;
        for (auto n : triangleSet) {
            if (trianglesAreNeighbors(triangles.back(), n)) {
                neighbors++;
                if (neighborsFull(n)) {
                    triangleSet.erase(n);
                }
            }

        }
        if (neighbors < targetNumberOfNeighbors(triangles.back())) {
            triangleSet.insert(triangles.back());
        }


        //        plot(c.first, k.first, st.first);

    };

    for (int i = 2; i < sorted.size(); i++) {
        auto c = sorted[i];
        auto top = stack.back();
        if (top.second == c.second) {//top and c in same chain
            auto k = top;
            stack.pop_back();
            while (!stack.empty() and stack.back().first->isRight(c.first, k.first) == c.second) {
                if (c.second == 1) {
                    addToTriangles(c, k, stack.back());
                } else {
                    addToTriangles(k, c, stack.back());
                }
                k = stack.back();
                stack.pop_back();
            }
            stack.push_back(k);
            stack.push_back(c);
        } else {
            while (!stack.empty()) {
                auto k = stack.back();
                stack.pop_back();
                if (!stack.empty()) {
                    if (k.second == 1) {

                        addToTriangles(c, k, stack.back());
                    } else {

                        addToTriangles(k, c, stack.back());

                    }
                }
            }
            stack.push_back(top);
            stack.push_back(c);

        }
    }
    return triangles;
}

auto inline connectAccrossDiagonals(std::vector<Edge> diagonals, int polygonId) {
    for (Edge d : diagonals) {

        std::vector<const Triangle*>*triangles1 = d.first->getTriangles(polygonId);
        std::vector<const Triangle*>*triangles2 = d.second->getTriangles(polygonId);

        std::vector<const Triangle*> opposing;
        for (auto it = triangles1->begin(); it != triangles1->end(); it++) {
            if (std::find(triangles2->begin(), triangles2->end(), *it) != triangles2->end()) {
                opposing.push_back(*it);
            }
        }

        for (int i = 0; i < 3; i++) {
            if (opposing[0]->points[i] != d.first and opposing[0]->points[i] != d.second) {

                opposing[0]->neighbours[i] = opposing[1];
                break;
            }
        }
        for (int i = 0; i < 3; i++) {
            if (opposing[1]->points[i] != d.first and opposing[1]->points[i] != d.second) {
                opposing[1]->neighbours[i] = opposing[0];
                break;
            }
        }
    }
}

bool inline isCW(std::vector<const Coords*> ring) {
    double sum = 0.0;
    for (int i = 0; i < ring.size(); i++) {
        sum += (ring[(i + 1) % ring.size()]->getX() - ring[i]->getX())*(ring[(i + 1) % ring.size()]->getY() + ring[i]->getY());
    }
    return sum > 0;
}

auto inline triangulatePolygon(std::vector<std::vector <const Coords*>> polygon, int polygonId) {
    std::vector<Triangle*> triangles;
    std::vector<Edge> diagonals{};
    auto monotones = splitToMonotone(polygon, &diagonals);
    //std::cout << "Monotones done" << std::endl;
    //plotPolygon(polygon, "black");
    for (auto mono : monotones) {
        /*if (isCW(mono)) {
            plotPolygon(std::vector<std::vector<const Coords*>>
            {
                mono
            }, "blue");

        } else {
            plotPolygon(std::vector<std::vector<const Coords*>>
            {
                mono
            }, "red");

        }*/
        auto t = triangulateMonotone(mono, polygonId);
        triangles.insert(triangles.end(), t.begin(), t.end());
    }
    //matplotlibcpp::show();

    connectAccrossDiagonals(diagonals, polygonId);
    return triangles;
}


#endif /* SRC_GEOMFUNC_H_ */
