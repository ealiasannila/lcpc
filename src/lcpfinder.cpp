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
#include<stack>

/*
 * Finds opposing by looking up intersection of the immidiate neighbours of either end of the base.
 * This can be either 1 or 2 vertices, of whcih 1 or 0 may be relevant. (irrelevant vertice is already inside funnel)
 * LOOKING UP INTERSECTION REQUIRES SORTING IMMIDIATE NEIGHBOURS. HOW TO MAKE FASTER:
 * ?-> have 
 * sorted vector of immidiate neighbours available. (increases insertion time to log N, now constant)
 */






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
        printTriangle2(*it);
        int ci = -1;
        for (int i = 0; i < 3; i++) {
            const Coords* neighbour = (*it)->points[i];
            if (neighbour != c) {
                insertToNset(nset, neighbour, fric, c);
            }
        }
        funnelQue.push_back(Funnel{c, *it});

    }
    return funnelQue;

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

void LcpFinder::findNeighboursInPolygon(const Coords* c, int polygon, nSet* nset) {

    std::deque<Funnel> funnelQue = this->initFQue(c, polygon, nset);
    std::pair<const Coords*, const Coords*> base;
    const Coords* o;
    Funnel f = funnelQue.front();
    while (!funnelQue.empty()) {
        f = funnelQue.front();
        //printFunnel(f);
        funnelQue.pop_front();
        f.stepForward(&funnelQue, nset, this->frictions[polygon], this->maxDist);
    }

    if (this->fences.find(polygon) != this->fences.end()) {
        for (auto p = nset->begin(); p != nset->end(); p++) {
            double fencecost = this->checkFences(polygon, c, p->first);
            if (c->flag == 2) {
                //std::cout << c->toString() << p->first->toString() << " cost increased by " << fencecost << std::endl;
            }
            p->second += fencecost;
        }
    }

}

nSet LcpFinder::findNeighbours(const Coords* c) {

    //std::cout<<"C: "<<c->toString()<<std::endl;
    nSet nset{};
    if (c != this->startPoint and c->flag == 1) {
        return nset;
    }
    std::vector<int> polygons = c->belongsToPolygons();
    for (int p : polygons) {
        if (p < 0) {
            std::cout << "HALOO!! -1 POLYGON\n";
        } else {
            findNeighboursInPolygon(c, p, &nset);
        }
    }

    if (c->flag == 2) {
        for (std::tuple<const Coords*, double, double> n : c->linearNeighbours) {
            //std::cout<<"linear fric: "<<n.second<<std::endl;
            insertToNset(&nset, std::get<0>(n), std::get<1>(n), c, 0);
        }
    }

    return nset;
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
    for(int i = 0; i<this->polygons.size();i++){
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
        if (node->flag == 1) {
            targetsFound++;
            if (targetsFound == this->numOfTargets) {
                std::cout << "break\n";
                break;
            }
        }


        minheap.pop();
        nSet neighbours = findNeighbours(node);

        for (std::pair<const Coords*, double> p : neighbours) {
            const Coords* n = p.first;
            if (n == 0) {
                std::cout << "neighbour is 0\n";
                exit(2);
            }
            double d = p.second;

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

int inSector(const Coords* a, const Coords* b, const Coords* right, const Coords* left) {
    int lob = b->isRight(a, left);
    int rob = b->isRight(a, right);

    if (lob == 1 and rob == -1) {
        return 1;
    }
    return -1;
}

bool CoordsInLinearNeighbour(const Coords* c, std::vector<std::tuple<const Coords*, double, double>>&ln) {
    for (std::tuple<const Coords *, double, double> t : ln) {
        if (std::get<0>(t) == c) {
            return true;
        }
    }
    return false;
}

const Coords* findEntry(const Coords* exit) {
    while (CoordsInLinearNeighbour(exit->getPred(), exit->linearNeighbours)) {
        exit = exit->getPred();
    }
    return exit;
}

bool checkCrossingAtPoint(const Coords* a, const Coords* b) {
    if (a->getPred() == 0) {
        return false;
    }

    const Coords* exitN1 = std::get<0>(a->linearNeighbours[1]);
    const Coords* exitN0 = std::get<0>(a->linearNeighbours[0]);

    int orientExit = exitN1->isRight(exitN0, a);
    int inSectorExit;
    int inSectorEntry;

    const Coords* entry = findEntry(a);

    const Coords* entryN1 = std::get<0>(entry->linearNeighbours[1]);
    const Coords* entryN0 = std::get<0>(entry->linearNeighbours[0]);
    int orientEntry = entryN1->isRight(entryN0, entry);

    if (orientExit == -1) {
        inSectorExit = inSector(a, b, exitN1, exitN0);
    } else if (orientExit == 0) {
        inSectorExit = b->isRight(exitN0, exitN1);
    } else {
        inSectorExit = inSector(a, b, exitN0, exitN1);
    }

    if (orientEntry == -1) {
        inSectorEntry = inSector(entry, entry->getPred(), entryN1, entryN0);
    } else if (orientEntry == 0) {
        inSectorEntry = entry->getPred()->isRight(entryN0, entryN1);
    } else {
        inSectorEntry = inSector(entry, entry->getPred(), entryN0, entryN1);
    }

    if (orientExit != 0) {
        inSectorExit *= orientExit;
    }
    if (orientEntry != 0) {
        inSectorEntry *= orientEntry;
    }

    if (true) {
        std::cout << std::endl;
        std::cout << "entry: " << entry->toString();
        std::cout << "exit: " << a->toString();
        std::cout << "B: " << b->toString();
        std::cout << "Pred: " << entry->getPred()->toString();
        std::cout << "exitN0" << exitN0->toString();
        std::cout << "exitN1" << exitN1->toString();
        std::cout << "entryN0" << entryN0->toString();
        std::cout << "entryN1" << entryN1->toString();
        std::cout << "orientExit: " << orientExit << std::endl;
        std::cout << "orientEntry: " << orientEntry << std::endl;
        std::cout << "inSectorExit: " << inSectorExit << std::endl;
        std::cout << "inSectorEntry: " << inSectorEntry << std::endl;

    }

    return inSectorExit != inSectorEntry;
}

double LcpFinder::checkFences(int polygon, const Coords* a, const Coords* b) {

    try {
        double crossingCost = 0;
        for (const Coords* fencePoint : this->fences.at(polygon)) {
            std::stack<const Coords*> stack;
            std::set<const Coords*> used;
            stack.push(fencePoint);
            while (!stack.empty()) {
                const Coords* fp = stack.top();
                used.insert(fp);
                stack.pop();
                for (std::tuple<const Coords*, double, double> p : fp->linearNeighbours) {

                    const Coords* next = std::get<0>(p);
                    if (used.find(next) != used.end()) {
                        continue;
                    }
                    int sc = segmentCrossing(fp, next, a, b);
                    if (sc == 1) {
                        //std::cout<<"  SC == 1"<<std::endl;
                        crossingCost += std::get<2>(p);
                    } else if (sc == 2 and a != fp and a != next and b != fp and b != next) {
                        crossingCost += std::get<2>(p);
                    }

                    /*
                            
                    if (a == fp and a->linearNeighbours.size() > 1) {
                        if (checkCrossingAtPoint(a, b)) {
                            //std::cout<<"  normal check cross"<<std::endl;
                            crossingCost += std::get<2>(p);
                        }
                    }
                     */
                    if (a == fp and !CoordsInLinearNeighbour(b, a->linearNeighbours)) {
                        crossingCost += std::get<2>(p) / 2; //LEAVING COST    


                    } else if (b == fp and !CoordsInLinearNeighbour(b, a->linearNeighbours)) {
                        crossingCost += std::get<2>(p) / 2;
                        //JOINING COST
                    }
                    /*if (next->belongsToPolygons().size() > 1) {
                        if (a == next and a->linearNeighbours.size() > 1 and checkCrossingAtPoint(a, b)) {
                            //std::cout<<"  last check cross";
                            crossingCost += std::get<2>(p);
                        }
                    }
                     */
                    if (next->belongsToPolygon(polygon)) {
                        stack.push(next);
                    }
                }
            }
        }
        return crossingCost;
    } catch (std::out_of_range & oor) {

        return 0;
    }
}

void LcpFinder::checkTargets(int polygon, Triangle * newTri) {
    try {
        std::forward_list<const Coords*>::const_iterator prev;
        for (auto it = this->targetPoints.at(polygon).cbegin(); it != this->targetPoints.at(polygon).cend();) {
            const Coords* target = *it;
            if (coordsInTriangle(newTri, target)) {
                newTri->interiorPoints.push_back(target);
                /*
                if (it == this->targetPoints.at(polygon).cbegin()) {
                    this->targetPoints.at(polygon).pop_front();
                    it = this->targetPoints.at(polygon).cbegin();
                } else {
                    it = this->targetPoints.at(polygon).erase_after(prev);
                }
                this->targetPoints[-1].push_front(target);
                 * */
                prev = it;
                it++;
            } else {
                prev = it;
                ++it;
            }
        }
    } catch (const std::out_of_range& oor) {
    }
}

void LcpFinder::subTriangles(Triangle* mainTriangle, int polygon, const Coords* c) {
    for (int i = 0; i < 3; i++) {
        Triangle* tri = new Triangle{};
        tri->points[0] = c;
        tri->points[1] = mainTriangle->points[i];
        tri->points[2] = mainTriangle->points[(i + 1) % 3];
        tri->neighbours[0] = mainTriangle->neighbours[(i + 2) % 3];
        std::cout << c->toString();
        std::cout << polygon << std::endl;
        c->addTriangle(tri, polygon);
        for (const Coords* interior : mainTriangle->interiorPoints) {
            if (coordsInTriangle(tri, interior)) {
                tri->interiorPoints.push_back(interior);
            }
        }
        this->checkTargets(polygon, tri);
        this->checkLinear(polygon, tri, false);
    }
}

void LcpFinder::checkLinear(int polygon, Triangle * newTri, bool subtri) {
    try {
        std::forward_list<const Coords*>::const_iterator prev;
        for (auto it = this->linePoints.at(polygon).cbegin(); it != this->linePoints.at(polygon).cend();) {
            const Coords* lp = *it;
            if (coordsInTriangle(newTri, lp)) {
                newTri->interiorPoints.push_back(lp);
                /*
                if (it == this->linePoints.at(polygon).cbegin()) {
                    this->linePoints.at(polygon).pop_front();
                    it = this->linePoints.at(polygon).cbegin();
                } else {
                    it = this->linePoints.at(polygon).erase_after(prev);
                }
                this->linePoints[-1].push_front(lp);
                this->subTriangles(newTri, polygon, lp);
                 * */
                if (subtri) {
                    this->subTriangles(newTri, polygon, lp);
                }
                prev = it;
                it++;

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

    std::vector<std::vector < p2t::Point*>> points(this->polygons.at(polygon).size());
    try {
        for (int ring = 0; ring<this->polygons.at(polygon).size(); ring++) {
            points[ring] = std::vector<p2t::Point*>(this->polygons.at(polygon).at(ring).size());
            for (int i = 0; i < points[ring].size(); i++) {
                const Coords* c = this->polygons.at(polygon).at(ring).at(i);
                points[ring][i] = new p2t::Point{c->getX(), c->getY()};
            }
        }
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
        std::cout << "triangle:\n";
        for (unsigned i = 0; i < 3; i++) {
            p2t::Point* point = triangle->GetPoint(i);
            c[i] = Coords(point->x, point->y);
            std::cout << c[i].toString();
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
            }
        }
        for (int n = 0; n < 3; n++) {
            if (triangles[i]->constrained_edge[n]) {
                const Coords* lb = newTri->points[(n + 2) % 3];
                const Coords* rb = newTri->points[(n + 1) % 3];
                double d = eucDistance(lb, rb);
                if (d>this->maxDist) {
                    std::cout << "intermediate points between" << lb->toString() << rb->toString();
                    std::cout << newTri;
                    int toAdd = std::ceil(d / this->maxDist) - 1;
                    std::cout << "adding: " << toAdd << std::endl;
                    for (int i = 1; i <= toAdd; i++) {
                        double x{((rb->getX() - lb->getX()) / (toAdd + 1)) * i + lb->getX()};
                        double y{((rb->getY() - lb->getY()) / (toAdd + 1)) * i + lb->getY()};
                        std::pair < std::tr1::unordered_set<Coords, CoordsHasher>::iterator, bool> success = this->coordmap.insert(Coords{x, y, polygon, false});

                        const Coords* intermediate = &*success.first;
                        intermediate->addToPolygon(polygon);
                        std::cout << intermediate->toString();

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

                        for (const Coords* interior : newTri->interiorPoints) {
                            if (coordsInTriangle(t1, interior)) {
                                t1->interiorPoints.push_back(interior);
                            }
                            if (coordsInTriangle(t2, interior)) {
                                t2->interiorPoints.push_back(interior);
                            }
                        }
                        checkTargets(polygon, t1); //change to using old triangle interior points to make faster
                        checkTargets(polygon, t2);

                        intermediate->addTriangle(t1, polygon);
                        intermediate->addTriangle(t2, polygon);
                    }
                }
            }
        }

        if (coordsInTriangle(newTri, this->startPoint)) {
            subTriangles(newTri, polygon, this->startPoint);
        }
        this->checkTargets(polygon, newTri);
        this->checkLinear(polygon, newTri, true);

    }
    for (std::vector<p2t::Point*> r : points) {
        for (p2t::Point* p : r) {

            delete p;
        }
    }

}

LcpFinder::~LcpFinder() {


}

void LcpFinder::addPolygon(std::vector<std::vector < p2t::Point*>> points, double friction) {
    //std::cout<<"adding POLYGON\n";
    int polygon = polygons.size();
    if (friction < this->minFriction) {
        minFriction = friction;
    }

    this->frictions.push_back(friction);
    this->triangulated.push_back(false);
    this->polygons.push_back(std::vector<std::vector<const Coords*>>(points.size()));

    for (int ringIndex = 0; ringIndex < points.size(); ringIndex++) {
        std::vector<p2t::Point*> ring = points[ringIndex];
        this->polygons.back()[ringIndex].reserve(ring.size());
        for (int i = 0; i < ring.size(); i++) {
            p2t::Point* point;
            if (ringIndex == 0) {
                point = ring[i];
            } else {
                point = ring[ring.size() - 1 - i];
            }
            std::pair < std::tr1::unordered_set<Coords>::iterator, bool> p = coordmap.insert(Coords(point->x, point->y, polygon, false));
            if (!p.second) {

                p.first->addToPolygon(polygon);
            }
            this->polygons.back()[ringIndex].push_back(&*p.first);
        }
    }
}

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

    std::pair < std::tr1::unordered_set<Coords>::iterator, bool> success = this->coordmap.insert(Coords(point->x, point->y, polygon, 2));
    if (!success.second) {
        success.first->addToPolygon(polygon);
    }
    const Coords* c = &*success.first;
    if (this->linePoints.find(polygon) == this->linePoints.end()) {
        this->linePoints[polygon] = std::forward_list<const Coords*>{c};
    } else {

        this->linePoints[polygon].push_front(c);
    }
    return c;
}

const Coords * LcpFinder::addLinePoint(std::array<double, 2>xy, int polygon) {

    std::pair < std::tr1::unordered_set<Coords>::iterator, bool> success = this->coordmap.insert(Coords(xy[0], xy[1], polygon, 2));
    if (!success.second) {
        success.first->addToPolygon(polygon);
    }
    const Coords* c = &*success.first;
    if (this->linePoints.find(polygon) == this->linePoints.end()) {
        this->linePoints[polygon] = std::forward_list<const Coords*>{c};
    } else {

        this->linePoints[polygon].push_front(c);
    }
    return c;
}

void LcpFinder::addBuffers(double d) {
    for (std::pair<int, std::vector<const Coords*>> polygon : this->fences) {
        for (const Coords* start : polygon.second) {
            std::stack<const Coords*> stack;
            std::set<const Coords*> used;
            stack.push(start);
            while (!stack.empty()) {
                const Coords* fp = stack.top();
                used.insert(fp);
                stack.pop();
                for (int i = 0; i < fp->linearNeighbours.size(); i++) {
                    std::tuple<const Coords*, double, double> n = fp->linearNeighbours[(i + 1) % fp->linearNeighbours.size()];
                    std::tuple<const Coords*, double, double> p = fp->linearNeighbours[i];
                    const Coords* next = std::get<0>(n);
                    const Coords* prev = std::get<0>(p);

                    addBufferPoint(prev, fp, next, d, polygon.first);

                    if (used.find(next) != used.end()) {
                        continue;
                    }
                    if (next->belongsToPolygons().size() == 1 and next->linearNeighbours.size() > 1) {
                        stack.push(next);
                    }
                }
            }
        }
    }
}

void LcpFinder::addBufferPoint(const Coords* prev, const Coords* c, const Coords* next, double d, int polygon) {

    if (prev == 0) {
        double t = d / eucDistance(c, next);
        double x = (1 - t) * c->getX() + t * next->getX();
        double y = (1 - t) * c->getY() + t * next->getY();
        this->addLinePoint(std::array<double, 2>{x, y}, polygon);
    }
    if (next == 0) {
        double t = d / eucDistance(prev, c)*-1;
        double x = (1 - t) * prev->getX() + t * c->getX();
        double y = (1 - t) * prev->getY() + t * c->getY();
        this->addLinePoint(std::array<double, 2>{x, y}, polygon);
    }
    double x2 = c->getX() - prev->getX();
    double y2 = c->getY() - prev->getY();
    double x1 = c->getX() - next->getX();
    double y1 = c->getY() - next->getY();

    double m = d / eucDistance(std::array<double, 2>{0, 0}, std::array<double, 2>{x1 + x2, y1 + y2});
    double lx = (x1 + x2) * m + c->getX();
    double ly = (y1 + y2) * m + c->getY();




    this->addLinePoint(std::array<double, 2>{lx, ly}, polygon);
    //this->addLinePoint(std::array<double, 2>{rx, ry}, polygon);



}

void LcpFinder::addLine(std::vector<p2t::Point*>* points, double friction, double crossing) {
    const Coords* prev = 0;

    for (auto it = points->begin(); it != points->end(); it++) {
        p2t::Point* p = *it;
        std::array<int, 2> polygons = containingPolygon(p);

        if (polygons[0] == -1) {
            continue;
        }
        const Coords* c = addLinePoint(p, polygons[0]);
        double minFric = friction; //std::min(this->frictions[polygons[0]], friction);
        if (polygons[1] != -1) {
            //minFric = std::min(this->frictions[polygons[1]], minFric);
            addLinePoint(p, polygons[1]);
        }
        if (prev != 0) {
            std::cout << c->toString();
            std::cout << prev->toString();
            int next;
            std::vector<const Coords*> crossings = segmentPolygonIntersection(polygons[0], prev, c, &next);
            for (int i = 0; i < crossings.size(); i++) {
                if (i == 0) {
                    prev->addLinearNeighbour(crossings[0], minFric, crossing);
                    crossings.front()->addLinearNeighbour(prev, minFric, crossing);
                    this->fences[polygons[0]].push_back(crossings[0]);
                }
                if (i == crossings.size() - 1) {
                    crossings.back()->addLinearNeighbour(c, minFric, crossing);
                    c->addLinearNeighbour(crossings.back(), minFric, crossing);

                }//ADD LINKS BETWEEN CROSSING POINTS CHOOSE MIN OF ALONG FRICTION AND THAT OF POLYGON
            }
            if (crossings.size() == 0) {
                prev->addLinearNeighbour(c, minFric, crossing);
                c->addLinearNeighbour(prev, minFric, crossing);
            }

            if (polygons.size() > 1 and crossing > 0) {
                this->fences[polygons[0]].push_back(c); //may not be at zero place
            }
        } else if (crossing > 0) {
            this->fences[polygons[0]].push_back(c);
        }
        prev = c;
    }
}

std::vector<const Coords*> LcpFinder::segmentPolygonIntersection(int polygon, const Coords* a, const Coords* b, int* nextPolygon) {
    std::vector<const Coords*> v;
    std::vector<std::vector < const Coords*>> p = this->polygons[polygon];
    for (int ring = 0; ring < p.size(); ring++) {
        for (int i = 0; i < p[ring].size(); i++) {
            const Coords* prev = p[ring][i];
            const Coords* next = p[ring][(i + 1) % p[ring].size()];
            if (segmentCrossing(a, b, prev, next) == 1) {
                std::array<double, 2> xy = lineIntersection(a, b, prev, next);
                std::pair < std::tr1::unordered_set<Coords>::iterator, bool> success = this->coordmap.insert(Coords{xy[0], xy[1], polygon});
                if (success.second) {

                    *nextPolygon = neighbouringPolygon(prev, next, polygon);
                    success.first->addToPolygon(*nextPolygon);
                }
                this->linePoints[polygon].push_front(&*success.first);
                this->linePoints[*nextPolygon].push_front(&*success.first);
                v.push_back(&*success.first);
            }
        }
    }
    return v;

}

/*
 * Returns at most two polygons inside which a point is.
 */

std::array<int, 2> LcpFinder::containingPolygon(p2t::Point * p) {
    std::array<int, 2>res{-1, -1};
    int found = 0;
    for (int i = 0; i<this->polygons.size(); i++) {
        std::vector<std::vector < const Coords*>> polygon = this->polygons[i];
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
    std::pair < std::tr1::unordered_set<Coords>::iterator, bool> success = this->coordmap.insert(Coords(startPoint->x, startPoint->y, polygon, 3));
    this->startPoint = & * success.first;
}
