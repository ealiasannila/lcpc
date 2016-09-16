/*
 * const Coords.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#include "Coords.h"
#include <sstream>
#include <math.h>
#include <iostream>
#include <algorithm>

Coords::Coords(double newx, double newy) {
    x = newx;
    y = newy;
    toStart = -1;
    this->predecessor = 0;
}

Coords::Coords(double newx, double newy, int polygon) {
    x = newx;
    y = newy;
    toStart = -1;
    this->addToPolygon(polygon);
    this->predecessor = 0;
}

Coords::Coords() {
    toStart = -1;
    x = -1;
    y = -1;
    this->predecessor = 0;
}

nContainer* Coords::getRightNeighbours(int polygon) const {
    try {
        return &this->rightNeighbours.at(polygon);
    } catch (const std::out_of_range& oor) {
        std::cout << "NoNeighbours: \n";
        std::cout << "polygon: " << polygon;
        exit(1);
    }
}

nContainer* Coords::getLeftNeighbours(int polygon) const {
    try {
        return &this->leftNeighbours.at(polygon);
    } catch (const std::out_of_range& oor) {
        std::cout << "NoNeighbours: \n";
        std::cout << "polygon: " << polygon;
        exit(1);
    }
}

allNeighIter Coords::getAllLeftN() const {
    return allNeighIter(this->leftNeighbours.begin(), this->leftNeighbours.end());
}

allNeighIter Coords::getAllRightN() const {
    return allNeighIter(this->rightNeighbours.begin(), this->rightNeighbours.end());
}

void Coords::addToPolygon(int polygon) const {
    this->leftNeighbours.insert(std::pair<int, nContainer>(polygon, nContainer()));
    this->rightNeighbours.insert(std::pair<int, nContainer>(polygon, nContainer()));
}

void Coords::addNeighbours(const Coords* l, const Coords* r, int polygon) const {
    try {


        /*
        auto ri = std::lower_bound(rightNeighbours.at(polygon).begin(), rightNeighbours.at(polygon).end(), r);
        rightNeighbours.at(polygon).insert(ri, r);
        
         */
        leftNeighbours.at(polygon).push_back(l);
        rightNeighbours.at(polygon).push_back(r);

        std::sort(leftNeighbours.at(polygon).begin(), leftNeighbours.at(polygon).end(), std::less<const Coords*>()); //... could this be avoided? Or is it that bad? they are small vectors...
//        std::sort(rightNeighbours.at(polygon).begin(), rightNeighbours.at(polygon).end(), std::less<const Coords*>()); //...



    } catch (const std::out_of_range& oor) {
        std::cout << "NoNeighbours: \n";
        std::cout << this->toString();
        std::cout << "polygon: " << polygon;
        exit(1);
    }
}

void Coords::sortNeighbours(unsigned int polygon) const {
    std::sort(leftNeighbours.at(polygon).begin(), leftNeighbours.at(polygon).end()); //... could this be avoided? Or is it that bad? they are small vectors...
    std::sort(rightNeighbours.at(polygon).begin(), rightNeighbours.at(polygon).end()); //...

}

void Coords::setToStart(double cost) const {
    toStart = cost;
}

void Coords::setPred(const Coords* pred) const {
    predecessor = pred;
}

std::string Coords::toString() const {
    std::stringstream sstm;
    sstm << "X: " << x << " Y: " << y << " Polygons: ";

    for (std::pair<int, nContainer> p : this->leftNeighbours) {
        sstm << p.first << " ";
    }
    return sstm.str();
}

int Coords::isRight(const Coords* c1, const Coords* c2) const {
    double d = (c2->getY() - c1->getY()) * (x - c2->getX()) - (c2->getX() - c1->getX()) * (y - c2->getY());
    if (d < 0) {
        return -1;
    }
    if (d > 0) {
        return 1;
    }
    return 0;
}

bool Coords::operator==(const Coords& c) const {
    double tol = 0.000001;
    if (std::abs(c.getX() - this->getX()) < tol and std::abs(c.getY() - this->getY()) < tol) {
        return true;
    }
    return false;
}

bool Coords::operator!=(const Coords& c) const {
    return !this->operator==(c);
}

Coords::~Coords() {
    // TODO Auto-generated destructor stub
}

