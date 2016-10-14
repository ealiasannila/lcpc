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


/*
 * Finds opposing by looking up intersection of the immidiate neighbours of either end of the base.
 * This can be either 1 or 2 vertices, of whcih 1 or 0 may be relevant. (irrelevant vertice is already inside funnel)
 * LOOKING UP INTERSECTION REQUIRES SORTING IMMIDIATE NEIGHBOURS. HOW TO MAKE FASTER:
 * ?-> have sorted vector of immidiate neighbours available. (increases insertion time to log N, now constant)
 */


const Coords* LcpFinder::getOpposing(const Coords* l, const Coords* r, int polygon) {
    nContainer intrsct;


    intrsct = intersection(l->getLeftNeighbours(polygon), r->getNeighbours(polygon));
    for (nContainer::iterator it = intrsct.begin(); it != intrsct.end(); it++) {
        if (*it != 0 and (*it)->isRight(l, r) == -1) {
            return *it;
        }
    }

    return 0;
}

std::deque<Funnel> LcpFinder::initFQue(const Coords* c, int polygon, nSet*neighbours) {
    nContainer* ln = c->getLeftNeighbours(polygon);
    if (ln->empty()) {
        this->triangulate(polygon);
    }

    nContainer* rn = c->getRightNeighbours(polygon);


    std::deque<Funnel> funnelQue{};
    nContainer::iterator rit = rn->begin();
    for (nContainer::iterator lit = ln->begin(); lit != ln->end(); lit++) {
        Funnel f(*lit, c, *rit);

        neighbours->insert(std::pair<const Coords*, int>(*lit, polygon));
        neighbours->insert(std::pair<const Coords*, int>(*rit, polygon));
        funnelQue.push_back(f);
        rit++;
    }
    return funnelQue;

}

void LcpFinder::findNeighboursInPolygon(const Coords* c, int polygon, nSet* neighbours) {
    std::deque<Funnel> funnelQue = initFQue(c, polygon, neighbours);
    std::pair<const Coords*, const Coords*> base;
    const Coords* o;
    while (!funnelQue.empty()) {
        Funnel f = funnelQue.front();
        funnelQue.pop_front();
        base = f.getBase();
        o = getOpposing(base.first, base.second, polygon);
        if (o != 0) {
            std::clock_t rb = std::clock();
            f.reactToOpposite(o, &funnelQue, neighbours, polygon);
            funnelQue.push_back(f);
        }

    }
}

nSet LcpFinder::findNeighbours(const Coords* c) {
    nSet neighbours{};
    std::vector<int> polygons = c->belongsToPolygons();
    for (int p : polygons) {
        findNeighboursInPolygon(c, p, &neighbours);
    }
    return neighbours;
}

double LcpFinder::toClosestEnd(const Coords* c) {
    double min = std::numeric_limits<double>::max();
    for (std::pair<int, std::vector < p2t::Point*>> polygon : targetPoints) {
        for (p2t::Point* target : polygon.second) {
            double d = eucDistance(target, c);
            if (d < min) {
                min = d;
            }
        }
    }
    return min * this->minFriction;
}

 struct cmpr {
        bool (*f)(const Coords* , const Coords* );

        cmpr(){
            f = &compAstar;
        }
        
        cmpr(bool (*comparefunction)(const Coords* , const Coords* )) {
            f = comparefunction;
        }

        bool operator()(const Coords* x, const Coords* y) const {
            return f(x, y);
        }
    };

std::deque<const Coords*> LcpFinder::leastCostPath(int algorithm) {
    int targetsFound = 0;
    this->startPoint2->setToStart(0);

    bool (*compareFunction)(const Coords* , const Coords* );
    

    if (algorithm == 0) {
        compareFunction = &compDijkstra;
    } else if (algorithm == 1) {
        compareFunction = &compAstar;
    }
    cmpr comparator{compareFunction};

    MinHeap<const Coords*, cmpr> minheap(comparator);



    minheap.push(this->startPoint2);

    int handled = 0;
    int old = 0;
    int total = this->coordmap.size();

    while (!minheap.empty()) {

        int percentage = handled * 100 / total;
        if (percentage != old) {
            old = percentage;
            std::cout << "\rSearching..." << percentage << "% done";
            fflush(stdout);
        }
        handled++;
        const Coords* node = minheap.top();
        if (node->target) {
            targetsFound++;
            if (targetsFound == this->numOfTargets) {
                break;
            }
        }
        
        
        minheap.pop();
        nSet neighbours = findNeighbours(node); 
        for (std::pair<const Coords*, int> p : neighbours) {
            const Coords* n = p.first;
            double d{node->getToStart() + eucDistance(node, n) * this->frictions[p.second]};
            if (n->getToStart() < 0) { // node has not yet been inserted into minheap
                n->setToStart(d);
                n->setToEnd(this->toClosestEnd(n));
                n->setPred(node);
                minheap.push(n);
            } else if (n->getToStart() > d) {
                n->setToStart(d);
                n->setPred(node);
                minheap.update(); //reorders minheap after changing priority of n
            }
        }
    }

    std::deque<const Coords*> res;
    for (std::pair<int, std::vector < p2t::Point*>> endpoints : this->targetPoints) {
        for (p2t::Point* ep : endpoints.second) {
            if (ep->x == this->startPoint2->getX() and ep->y == this->startPoint2->getY()) {
                continue;
            }
            std::tr1::unordered_set<Coords>::iterator it = this->coordmap.find(Coords(ep->x, ep->y));
            if (it != this->coordmap.end()) {
                res.push_back(&*it);
            }
        }
    }
    return res;
}

void LcpFinder::triangulate(int polygon) {

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
    bool stp = true;
    std::vector<p2t::Point*> steinerPoints;
    try {
        steinerPoints = this->targetPoints.at(polygon);
    } catch (const std::out_of_range& oor) {
        stp = false;
    }
    if (stp) {
        for (p2t::Point* sp : steinerPoints) {
            cdt.AddPoint(sp);
        }
    }
    cdt.Triangulate();
    std::vector<p2t::Triangle*> triangles = cdt.GetTriangles();
    for (std::vector<p2t::Triangle*>::iterator it = triangles.begin(); it != triangles.end(); it++) {
        p2t::Triangle* triangle = *it;
        const Coords * cp[3];
        Coords c[3];
        for (unsigned i = 0; i < 3; i++) {
            p2t::Point* point = triangle->GetPoint(i);
            c[i] = Coords(point->x, point->y);
        }
        // if triangle orientation is clockwise turn it to CCW
        if (c[0].isRight(&c[1], &c[2]) == 1) {
            Coords tmp = c[2];
            c[2] = c[1];
            c[1] = tmp;
        }
        for (unsigned i = 0; i < 3; i++) {
            std::tr1::unordered_set<Coords>::iterator f = coordmap.find(c[i]);
            if (f != coordmap.end()) {
                cp[i] = &*f;
            } else {
            }
        }
        for (unsigned i = 0; i < 3; i++) {
            int l = i - 1;
            int r = i + 1;
            if (l == -1) {
                l = 2;
            }
            if (r == 3) {
                r = 0;
            }

            cp[i]->addNeighbours(cp[l], cp[r], polygon);

        }
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
    for (std::pair<int, std::vector < p2t::Point*>> polygon : targetPoints) {
        for (p2t::Point* point : polygon.second) {
            delete point;
        }
    }

}

void LcpFinder::addPolygon(std::vector<std::vector<p2t::Point*>> points, double friction) {
    if (friction < this->minFriction) {
        minFriction = friction;
    }
    this->frictions.push_back(friction);
    this->polygons.push_back(points);
    for (std::vector<p2t::Point*> ring : points) {
        for (p2t::Point* point : ring) {
            std::pair < std::tr1::unordered_set<Coords>::iterator, bool> p = coordmap.insert(Coords(point->x, point->y, polygons.size() - 1, false));
            if (!p.second) {
                p.first->addToPolygon(polygons.size() - 1);
            }

        }
    }
}

void LcpFinder::addSteinerPoint(p2t::Point* steinerpoint, int polygon) {
    if (this->targetPoints.find(polygon) == this->targetPoints.end()) {
        this->targetPoints[polygon] = std::vector<p2t::Point*>{steinerpoint};
    } else {
        this->targetPoints[polygon].push_back(steinerpoint);
    }
    std::pair < std::tr1::unordered_set<Coords>::iterator, bool> success = this->coordmap.insert(Coords(steinerpoint->x, steinerpoint->y, polygon, true));
    if (!success.second) {
        std::cout << "DID NOT INSERT\n";
        std::cout << success.first->toString() << std::endl;
        std::cout << Coords(steinerpoint->x, steinerpoint->y, polygon, true).toString() << std::endl;
    } else {
        this->numOfTargets++;
    }
}

void LcpFinder::addStartPoint(p2t::Point* startPoint, int polygon) {
    this->addSteinerPoint(startPoint, polygon);
    this->startPoint2 = & * this->coordmap.find(Coords(startPoint->x, startPoint->y));

}
