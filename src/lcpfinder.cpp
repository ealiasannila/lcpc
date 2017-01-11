/*
 * lcpfinder.cpp
 *
 *  Created on: Aug 25, 2016
 *      Author: elias
 */

#include "lcpfinder.h"
#include <algorithm>
#include <stdexcept>      // std::out_of_range
#include <iomanip>
#include <iostream>
#include<boost/heap/fibonacci_heap.hpp>

/*
 * Finds opposing by looking up intersection of the immidiate neighbours of either end of the base.
 * This can be either 1 or 2 vertices, of whcih 1 or 0 may be relevant. (irrelevant vertice is already inside funnel)
 * LOOKING UP INTERSECTION REQUIRES SORTING IMMIDIATE NEIGHBOURS. HOW TO MAKE FASTER:
 * ?-> have sorted vector of immidiate neighbours available. (increases insertion time to log N, now constant)
 */



void insertToNset(nSet* nset, std::pair<const Coords*, double> pair) {
    auto p = nset->insert(pair);
    if (!p.second and p.first->second > pair.second) {
        p.first->second = pair.second;
    }
}

void insertToNset(nSet* nset, const Coords* c, double fric) {
    insertToNset(nset, std::make_pair(c, fric));
}

void LcpFinder::setMaxD(double d) {
    this->maxDist = d;
}

std::deque<Funnel> LcpFinder::initFQue(const Coords* c, int polygon, nSet*nset) {
    if (c == 0) {
        std::cout << "WTF<!!" << std::endl;
        exit(1);
    }
    if (!this->triangulated[polygon]) {
        this->triangulate(polygon);
    }
    std::deque<Funnel> funnelQue{};
    double fric = this->frictions[polygon];
    for (auto it = c->getTriangles(polygon)->begin(); it != c->getTriangles(polygon)->end(); it++) {
        int ci = -1;
        for (int i = 0; i < 3; i++) {
            const Coords* neighbour = (*it)->points[i];
            if (neighbour != c) {
                insertToNset(nset, neighbour, fric);
            }
        }
        funnelQue.push_back(Funnel{c, *it});

    }
    return funnelQue;

}

void printTriangle2(const Triangle* t) {
    std::cout << std::fixed;
    std::cout << "Triangle : " << t << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "p: " << t->points[i]->getX() << "," << t->points[i]->getY() << std::endl;
        std::cout << "n: " << t->neighbours[i] << std::endl;
    }
}

void printFunnel(Funnel f) {
    std::cout << "f: " << f.apex;
    std::cout << "apex: " << f.getApex()->getX() << "," << f.getApex()->getY() << std::endl;

    std::cout << "right: " << f.firstRight->getX() << "," << f.firstRight->getY() << std::endl;
    std::cout << "left: " << f.firstLeft->getX() << "," << f.firstLeft->getY() << std::endl;
    std::cout << "base: " << f.base << std::endl;
    std::cout << "tri: ";
    printTriangle2(f.t);
    std::cout << std::endl;
    if (f.getOpposing() != 0) {
        std::cout << "opposing: " << f.getOpposing()->getX() << "," << f.getOpposing()->getY() << std::endl;
    } else {
        std::cout << "NO OPPOSING" << std::endl;
    }
}

void LcpFinder::findNeighboursInPolygon(const Coords* c, int polygon, nSet* neighbours) {
    std::deque<Funnel> funnelQue = this->initFQue(c, polygon, neighbours);

    std::pair<const Coords*, const Coords*> base;
    const Coords* o;
    Funnel f = funnelQue.front();
    while (!funnelQue.empty()) {
        f = funnelQue.front();
        //printFunnel(f);
        funnelQue.pop_front();
        f.stepForward(&funnelQue, neighbours, this->frictions[polygon], this->maxDist);
    }
}

nSet LcpFinder::findNeighbours(const Coords* c) {

    //std::cout<<"C: "<<c->toString()<<std::endl;
    nSet neighbours{};
    if (c != this->startPoint and c->target) {
        return neighbours;
    }
    std::vector<int> polygons = c->belongsToPolygons();
    for (int p : polygons) {
        if (p != -1) {
            findNeighboursInPolygon(c, p, &neighbours);
        }
    }

    if (c->linePt) {
        std::vector<std::pair <const Coords*, double>>*lineNeighbours = c->getNeighbours(-1);
        for (auto it = lineNeighbours->begin(); it != lineNeighbours->end(); it++) {
            insertToNset(&neighbours, *it);
        }
    }
    return neighbours;
}

double LcpFinder::toClosestEnd(const Coords* c) {
    double min = std::numeric_limits<double>::max();
    for (std::pair<int, std::forward_list < const Coords*>> polygon : targetPoints) {
        for (const Coords* target : polygon.second) {
            double d = eucDistance(target, c);
            if (d < min) {
                min = d;
            }
        }
    }
    return min * this->minFriction;
}

std::deque<const Coords*> LcpFinder::leastCostPath(int algorithm) {
    if (this->startPoint == 0) {
        std::cout << "No startpoint set in finder. Is the startpoint inside any polygon?\n" << std::endl;
        exit(1);
    }

    for (int i = 0; i<this->polygons.size(); i++) {
        this->triangulate(i);
    }

    int targetsFound = 0;
    this->startPoint->setToStart(0);
    /*
    bool (*compareFunction)(const Coords*, const Coords*);
    

    if (algorithm == 0) {
        compareFunction = &compDijkstra;
    } else if (algorithm == 1) {
        compareFunction = &compAstar;
    }
    cmpr comparator{compareFunction};
     */
    //MinHeap<const Coords*, cmpr> minheap(comparator);
    //MinHeap<const Coords*, cmpr> minheap(comparator);
    boost::heap::fibonacci_heap<const Coords*, boost::heap::compare < Coords::cmpr>> minheap;

    this->startPoint->handle = minheap.push(this->startPoint);

    int handled = 0;
    int old = 0;
    int total = this->coordmap.size();
    int percent = total / 100;

    while (!minheap.empty()) {

        if (handled - old > percent) {
            int percentage = handled * 100 / total;
            old = handled;
            std::cout << "\rSearching..." << percentage << "% done";
            fflush(stdout);
        }
        handled++;

        const Coords* node = minheap.top();
        //std::cout<<node->toString()<<std::endl;
        if (node == 0) {
            std::cout << "current node = 0\n";
            exit(2);
        }
        if (node->target) {
            targetsFound++;
            if (targetsFound == this->numOfTargets) {
                break;
            }
        }


        minheap.pop();
        nSet neighbours = findNeighbours(node);

        //std::cout << node->toString() << "# of neighbours: " << neighbours.size() << std::endl;
        for (std::pair<const Coords*, int> p : neighbours) {
            const Coords* n = p.first;
            if (n == 0) {
                std::cout << "neighbour is 0\n";
                exit(2);
            }
            double d{node->getToStart() + eucDistance(node, n) * p.second};
            if (n->getToStart() < 0) { // node has not yet been inserted into minheap
                n->setToStart(d);
                //n->setToEnd(this->toClosestEnd(n));
                n->setPred(node);
                n->handle = minheap.push(n);
            } else if (n->getToStart() > d) {
                n->setToStart(d);
                n->setPred(node);
                minheap.increase(n->handle);
                //minheap.decrease(n); //reorders minheap after changing priority of n
            }
        }

    }

    std::deque<const Coords*> res{};
    for (std::pair<int, std::forward_list < const Coords*>> polygon : this->targetPoints) {

        for (const Coords* ep : polygon.second) {
            if (ep == 0) {
                std::cout << "EP == 0\n";
                continue;
            }

            if (ep == this->startPoint) {
                continue;
            }
            res.push_back(ep);

        }
    }
    return res;
}

void LcpFinder::checkTargets(int polygon, Triangle * newTri) {
    try {
        std::forward_list<const Coords*>::const_iterator prev;
        for (auto it = this->targetPoints.at(polygon).cbegin(); it != this->targetPoints.at(polygon).cend();) {

            const Coords* target = *it;


            if (coordsInTriangle(newTri, target)) {
                newTri->interiorPoints.push_back(target);
                if (it == this->targetPoints.at(polygon).cbegin()) {
                    this->targetPoints.at(polygon).pop_front();
                    it = this->targetPoints.at(polygon).cbegin();
                } else {
                    it = this->targetPoints.at(polygon).erase_after(prev);
                }
                this->targetPoints[-1].push_front(target);
            } else {
                prev = it;
                ++it;
            }
        }
    } catch (const std::out_of_range& oor) {
    }
}

void LcpFinder::triangulate(int polygon) {
    if (this->triangulated[polygon]) {
        return;
    }
    this->triangulated[polygon] = true;

    std::vector<std::vector < p2t::Point*>> points;
    try {
        points = this->polygons.at(polygon);
    } catch (const std::out_of_range& oor) {
        std::cout << "EXCEPTION WHEN GETTING POLYGON\n";
        std::cout << "polygon:" << polygon << std::endl;
        exit(1);
    }
    p2t::CDT cdt{points[0]}; //Constrained Delaunay triangulator with outer edges.
    for (unsigned int hole = 1; hole < points.size(); hole++) {
        cdt.AddHole(points[hole]);
    }



    cdt.Triangulate();
    std::vector<p2t::Triangle*> triangles = cdt.GetTriangles();
    std::vector<Triangle*> newTriangles{};
    for (std::vector<p2t::Triangle*>::iterator it = triangles.begin(); it != triangles.end(); it++) {
        p2t::Triangle* triangle = *it;
        std::array<const Coords*, 3> cp;
        Coords c[3];
        for (unsigned i = 0; i < 3; i++) {
            p2t::Point* point = triangle->GetPoint(i);
            c[i] = Coords(point->x, point->y);
        }
        // if triangle orientation is clockwise turn it to CCW
        if (c[0].isRight(&c[1], &c[2]) == 1) {
            std::cout << "CLOCKWISE\n";
            Coords tmp = c[2];
            c[2] = c[1];
            c[1] = tmp;
        }

        for (unsigned i = 0; i < 3; i++) {
            std::tr1::unordered_set<Coords>::iterator f = coordmap.find(c[i]);
            if (f != coordmap.end()) {
                cp[i] = &*f;
            } else {
                std::cout << "Triangulated point doesn't match point in coordmap. Terminating\n";
                exit(1);
            }
        }
        Triangle* t = new Triangle{};
        t->points = cp;
        for (int i = 0; i < 3; i++) {
            cp[i]->addTriangle(t, polygon);
        }
        newTriangles.push_back(t);
    }
    for (int i = 0; i < triangles.size(); i++) {
        Triangle* newTri = newTriangles[i];
        for (int n = 0; n < 3; n++) {
            if (!triangles[i]->constrained_edge[n]) {
                auto it = std::find(triangles.begin(), triangles.end(), triangles[i]->GetNeighbor(n));
                if (it == triangles.end()) {
                    std::cout << "NEIGHBOUR NOT FOUND!!" << std::endl;
                } else {
                    int index = std::distance(triangles.begin(), it);
                    Triangle* newNeighbour = newTriangles.at(index);
                    newTri->neighbours[n] = newNeighbour;
                }
            } else {
                //DO THE INTERMEDIATE POINTS
                const Coords* lb = newTri->points[(n + 2) % 3];
                const Coords* rb = newTri->points[(n + 1) % 3];
                double d = eucDistance(lb, rb);
                if (d>this->maxDist) {
                    int toAdd = std::ceil(d / this->maxDist) - 1;
                    for (int i = 1; i <= toAdd; i++) {
                        double x{((rb->getX() - lb->getX()) / (toAdd + 1)) * i + lb->getX()};
                        double y{((rb->getY() - lb->getY()) / (toAdd + 1)) * i + lb->getY()};
                        std::pair<std::tr1::unordered_set<Coords, CoordsHasher>::iterator, bool> success =  this->coordmap.insert(Coords{x, y, polygon, false});
                        const Coords* intermediate = &*success.first;
                        intermediate->addToPolygon(polygon);
                        newTri->interiorPoints.push_back(intermediate);

                        const Coords* thirdCorner = newTri->points[n];
                        Triangle* t1 = new Triangle{};
                        t1->points[0] = intermediate;
                        t1->points[1] = thirdCorner;
                        t1->points[2] = rb;
                        Triangle* t2 = new Triangle{};
                        t2->points[0] = intermediate;
                        t2->points[1] = lb;
                        t2->points[2] = thirdCorner;

                        t1->neighbours[2] = t2;
                        t1->neighbours[0] = newTri->neighbours[(n + 2) % 3];

                        t2->neighbours[1] = t1;
                        t2->neighbours[0] = newTri->neighbours[(n + 1) % 3];
                        
                        intermediate->addTriangle(t1, polygon);
                        intermediate->addTriangle(t2, polygon);
                    }
                }
            }
        }
        if (coordsInTriangle(newTri, this->startPoint)) {
            for (int i = 0; i < 3; i++) {
                Triangle* tri = new Triangle{};
                tri->points[0] = this->startPoint;
                tri->points[1] = newTri->points[i];
                tri->points[2] = newTri->points[(i + 1) % 3];
                tri->neighbours[0] = newTri->neighbours[(i + 2) % 3];
                this->startPoint->addTriangle(tri, polygon);
                this->checkTargets(polygon, tri);
            }
        }
        this->checkTargets(polygon, newTri);

    }
}

LcpFinder::~LcpFinder() {
    for (std::vector<std::vector < p2t::Point*>> polygon : polygons) {
        for (std::vector<p2t::Point*> ring : polygon) {
            for (p2t::Point* point : ring) {
                delete point;
            }
        }
    }
    for (std::pair<int, std::forward_list <const Coords*>> polygon : targetPoints) {
        for (const Coords* point : polygon.second) {
            delete point;
        }
    }

}

void LcpFinder::addPolygon(std::vector<std::vector < p2t::Point*>> points, double friction) {
    //std::cout<<"adding POLYGON\n";
    int polygon = polygons.size();
    if (friction < this->minFriction) {
        minFriction = friction;
    }

    this->frictions.push_back(friction);
    this->triangulated.push_back(false);
    this->polygons.push_back(points);

    bool ext = true;

    for (std::vector<p2t::Point*> ring : points) {
        ext = false;
        const Coords* prev = 0;
        const Coords* first = 0;
        for (int i = 0; i < ring.size(); i++) {
            p2t::Point* point;
            if (ext) {
                point = ring[i];
            } else {
                point = ring[ring.size() - 1 - i];
            }
            std::pair < std::tr1::unordered_set<Coords>::iterator, bool> p = coordmap.insert(Coords(point->x, point->y, polygon, false));
            if (!p.second) {
                p.first->addToPolygon(polygon);
            }

            const Coords* cur = &*p.first;
            if (prev != 0) {
                prev->addNeighbours(cur, polygon, friction);
                cur->addNeighbours(prev, polygon, friction);
            } else {
                first = cur;
            }
            prev = &*p.first;

        }
        prev->addNeighbours(first, polygon, friction);
        first->addNeighbours(prev, polygon, friction, true);
    }
}

/*
void LcpFinder::addLine(std::vector<p2t::Point*>* line, double frictionForwards, double frictionBackwards, std::array<int, 2> polygons) {
    const Coords* prev = 0;
    for (int i = 0; i < line->size(); i++) {
        Coords f{line->at(i)->x, line->at(i)->y};
        auto it = this->coordmap.find(f);
        if (it == this->coordmap.end()) {
            std::cout << f.toString();
            const Coords* closest = 0;
            const Coords* c;
            for (auto it = this->getCoordmap()->begin(); it != this->getCoordmap()->end(); it++) {
                c = &*it;
                if (closest == 0 or eucDistance(c, &f) < eucDistance(closest, &f)) {
                    closest = c;
                }
            }
            std::cout << closest->toString();
            std::cout << "A linear feature's point was not found from cost surface!\n";
            exit(1);
        }
        const Coords* cur = &*it;
        for (int pol : polygons) {
            if (pol > 0 and pol<this->polygons.size()) {
                if (prev != 0) {
                    cur->addNeighbours(prev, pol, frictionBackwards);
                    prev->addNeighbours(cur, pol, frictionForwards);
                }
            }
        }


    }
}
 */
void LcpFinder::addTargetPoint(p2t::Point* targetpoint, int polygon) {
    std::pair < std::tr1::unordered_set<Coords>::iterator, bool> success = this->coordmap.insert(Coords(targetpoint->x, targetpoint->y, polygon, true));
    if (!success.second) {
        std::cout << "DID NOT INSERT steinerpoint: " << targetpoint->x << "," << targetpoint->y << "into polygon " << polygon << "\n";
        std::cout << success.first->toString() << std::endl;
        exit(1);
    } else {
        this->numOfTargets++;
    }
    const Coords* targetCoords = &*success.first;
    if (this->targetPoints.find(polygon) == this->targetPoints.end()) {
        this->targetPoints[polygon] = std::forward_list<const Coords*>{targetCoords};
    } else {
        this->targetPoints[polygon].push_front(targetCoords);
    }

}

const Coords * LcpFinder::addLinePoint(p2t::Point* point, int polygon) {

    if (this->linePoints.find(polygon) == this->linePoints.end()) {
        this->linePoints[polygon] = std::vector<p2t::Point*>{point};
    } else {
        this->linePoints[polygon].push_back(point);
    }
    std::pair < std::tr1::unordered_set<Coords>::iterator, bool> success = this->coordmap.insert(Coords(point->x, point->y, polygon, false, true));
    if (!success.second) {
        success.first->addToPolygon(polygon);
    }
    return &*success.first;
}

void LcpFinder::addLine(std::vector<p2t::Point*>* points, double friction) {
    const Coords* prev = 0;
    for (auto it = points->begin(); it != points->end(); it++) {
        p2t::Point* p = *it;
        std::array<int, 2> polygons = containingPolygon(p);
        const Coords* c = addLinePoint(p, polygons[0]);
        if (polygons[1] != -1) {
            addLinePoint(p, polygons[1]);
        }
        c->addToPolygon(-1);
        if (prev != 0) {
            c->addNeighbours(prev, -1, friction);
            prev->addNeighbours(c, -1, friction);

        }
        prev = c;
    }
}

/*
 * Returns at most two polygons inside which a point is.
 */

std::array<int, 2> LcpFinder::containingPolygon(p2t::Point * p) {
    std::array<int, 2>res{-1, -1};
    int found = 0;
    for (int i = 0; i<this->polygons.size(); i++) {
        std::vector<std::vector < p2t::Point*>> polygon = this->polygons[i];
        if (inside(polygon, p)) {
            res[found] = i;
            found++;
            if (found == 2) {
                break;
            }
        }
    }
    return res;
}

void LcpFinder::addStartPoint(p2t::Point* startPoint, int polygon) {
    std::pair < std::tr1::unordered_set<Coords>::iterator, bool> success = this->coordmap.insert(Coords(startPoint->x, startPoint->y, polygon, true));
    this->startPoint = & * success.first;

}
