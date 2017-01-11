/*
 * Funnel.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#include "funnel.h"
#include <sstream>
#include <iostream>
#include <algorithm>
#include "geomfunc.h"

Funnel::Funnel(const Coords* a, const Triangle* t) {
    for (int i = 0; i < 3; i++) {
        if (t->points[i] == a) {
            this->apex = a;
            this->base = i;
            this->firstLeft = t->points[(i + 2) % 3];
            this->firstRight = t->points[(i + 1) % 3];
            this->t = t;
            return;
        }
    }
}

Funnel::Funnel(const Coords* fl, const Coords* fr, const Coords* a, const Triangle* t, int b) {
    this->apex = a;
    this->base = b;
    this->firstLeft = fl;
    this->firstRight = fr;
    this->t = t;

}

//TODO varaudu jos deque on tyhjä tai jos törmää apexiin?

const Coords* Funnel::getOpposing() {
    const Triangle* nt = this->t->neighbours[this->base];
    if (nt == 0) {
        return 0;
    }
    for (int i = 0; i < 3; i++) {
        if (nt->points[i] == this->t->points[(this->base + 1) % 3]) {
            return nt->points[(i + 1) % 3];
        }
    }
    return 0;
}

std::string Funnel::toString() {
    std::stringstream sstm;
    sstm << "a: " << this->apex;
    return sstm.str();

}

/*
 * o is inside first segments of funnel, in which case funnel is split
 * so that left and right chains form two new funnels with the other chain
 * formed by apex and c.
 */
int Funnel::getNewBase(const Triangle* nt, const Coords* oldOpposing, int numOfSteps) {
    for (int i = 0; i < 3; i++) {
        if (nt->points[i] == oldOpposing) {
            return (i + numOfSteps) % 3;
        }
    }
}

Funnel Funnel::split() {
    const Coords* oldOpposing = this->getOpposing();
    const Coords* oldRight = this->firstRight;
    const Triangle* rightTri = this->t->neighbours[this->base];
    int rightBase;
    if (rightTri != 0) {
        rightBase = this->getNewBase(rightTri, oldOpposing, 1);
    }
    const Triangle* leftTri = this->t->neighbours[this->base];
    if (leftTri != 0) {
        this-> base = this->getNewBase(leftTri, oldOpposing, 2);
        this-> t = leftTri;
        this->firstRight = oldOpposing;
    }
    return Funnel(oldOpposing, oldRight, this->apex, rightTri, rightBase);
}

/*tests if o is inside first sector:
 * 0 -> in
 * -1 -> left
 * 1 -> right
 */
int Funnel::inFirstSector(const Coords* o) {
    int lo = o->isRight(this->apex, this->firstLeft);
    int ro = o->isRight(this->apex, this->firstRight);

    if (lo != 1) {
        return -1;
    }
    if (ro != -1) {
        return 1;
    }
    return 0;
}

void Funnel::intermediatesAtBase(double maxD, nSet* nset, double friction) {
    const Coords* baseL = this->t->points[(this->base + 1) % 3];
    const Coords* baseR = this->t->points[(this->base + 2) % 3];

    std::array<double, 2> intersectionLeft = lineIntersection(this->apex, this->firstLeft, baseL, baseR);
    std::array<double, 2> intersectionRight = lineIntersection(this->apex, this->firstRight, baseL, baseR);
    std::pair<int, const Triangle*> common = commonTriangle(baseL, baseR, this->t);

    //CHECK IF INTERSECTION IS AT EITHER END OF BASE!
    const Coords* epl = new Coords{intersectionLeft[0], intersectionLeft[1]};
    const Coords* epr = new Coords{intersectionRight[0], intersectionRight[1]};
    double d = eucDistance(epl, epr);
    int toAdd = std::ceil(d / maxD) - 1;
    for (int i = 1; i <= toAdd; i++) {
        double x{((epr->getX() - epl->getX()) / (toAdd + 1)) * i + epl->getX()};
        double y{((epr->getY() - epl->getY()) / (toAdd + 1)) * i + epl->getY()};

        //delete when popping from heap in main finder function.
        const Coords* inter = new Coords{x, y, common.first, false, 9};
        // DO FOLLOWING ALSO FOR EPR AND EPL UNLESS THEY ARE BASE SEE ABOVE
        int third = thirdCorner(common.second, baseL, baseR);
        const Coords* thirdCorner = common.second->points[third];
        inter->addToPolygon(common.first);


        Triangle* t1 = new Triangle{};
        t1->points[0] = inter;
        t1->points[1] = thirdCorner;
        t1->points[2] = baseL;
        Triangle* t2 = new Triangle{};
        t2->points[0] = inter;
        t2->points[1] = baseR;
        t2->points[2] = thirdCorner;

        t1->neighbours[2] = t2;
        t1->neighbours[0] = common.second->neighbours[(third + 2) % 3];

        t1->neighbours[1] = t2;
        t1->neighbours[0] = common.second->neighbours[(third + 2) % 3];

        inter->addTriangle(t1, common.first);
        inter->addTriangle(t2, common.first);
        std::pair<const Coords*, int> pair {inter, friction};
        nset->insert(pair);
        
    }

}

/*
 * Returns location of opposite node from the base of the funnel:
 * Does one of the following depending on the location of opposing vertex o:
 *  1. splits funnel
 *  2. shrinks either chain
 *  3. expands either chain
 */
void Funnel::stepForward(std::deque<Funnel>* funnelQueue, nSet* nset, double friction, double maxD) {

    for (const Coords* interior : this->t->interiorPoints) {
        if (inFirstSector(interior) == 0) {
            auto res = nset->insert(std::pair<const Coords*, int>(interior, friction));
            if (!res.second and res.first->second > friction) {
                res.first->second = friction;
            }
        }
    }
    const Coords* oldOpposing = this->getOpposing();
    if (oldOpposing == 0) {
        //this->intermediatesAtBase(maxD, nset, friction);
        return;
    }
    int inFirst = this->inFirstSector(oldOpposing);
    if (inFirst == 0) {
        Funnel n = this->split();
        if (n.t != 0) {
            funnelQueue->push_back(n);
        }
        if (this->t != 0) {
            funnelQueue->push_back(*this);
        }
        auto res = nset->insert(std::pair<const Coords*, int>(oldOpposing, friction));
        if (!res.second and res.first->second > friction) {
            res.first->second = friction;
        }
    } else {
        int numOfSteps;
        if (inFirst == -1) {
            numOfSteps = 1;
        } else {
            numOfSteps = 2;
        }
        this->t = this->t->neighbours[base];
        if (t != 0) {
            this->base = this->getNewBase(this->t, oldOpposing, numOfSteps);
            funnelQueue->push_back(*this);
        }
    }
}

Funnel::~Funnel() {
    // TODO Auto-generated destructor stub
}

