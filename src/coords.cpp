/*
 * const Coords.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#include "coords.h"
#include <sstream>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <limits>
#include <iomanip>
#include <functional>

Coords::Coords(double newx, double newy) {
    x = newx;
    y = newy;
    toStart = -1;
    this->predecessor = 0;
    this->flag = -1;
}

Coords::Coords(double newx, double newy, int polygon, int flag) {
    x = newx;
    y = newy;
    toStart = -1;
    this->addToPolygon(polygon);
    this->predecessor = 0;
    this->flag = flag;
}

Coords::Coords(double newx, double newy, int polygon) {
    x = newx;
    y = newy;
    toStart = -1;
    this->addToPolygon(polygon);
    this->predecessor = 0;
    this->flag = 0;
}

Coords::Coords() {
    toStart = -1;
    x = -1;
    y = -1;
    this->predecessor = 0;
    this->flag = -1;
}

std::vector<const Triangle*>*Coords::getTriangles(int polygon) const {
    return &(this->triangles.at(polygon));
}

void Coords::addTriangle(Triangle* t, int p) const {
    this->triangles.at(p).push_back(t);
}

std::vector<int> Coords::belongsToPolygons() const {
    std::vector<int> polygons;
    for (std::pair<int, std::vector<const Triangle*>> p : this->triangles) {

        polygons.push_back(p.first);
    }
    return polygons;
}

bool Coords::belongsToPolygon(int p) const {
    return this->triangles.find(p) != this->triangles.end();
}

void Coords::addToPolygon(int polygon) const {
    this->triangles.insert(std::make_pair(polygon, std::vector<const Triangle*>{}));
}

void Coords::addLinearNeighbour(const Coords* n, double friction) const{
    this->linearNeighbours.push_back(std::make_pair(n, friction));
}

void Coords::setToStart(double cost) const {
    this->toStart = cost;
}

void Coords::setToEnd(double cost) const {
    this->toEnd = cost;
}

void Coords::setPred(const Coords* pred) const {
    predecessor = pred;
}

std::string Coords::toString() const {
    std::stringstream sstm;
    sstm << std::setprecision(4);
    sstm << std::fixed;
    sstm << "xy: " << this->x << "," << this->y << std::endl;

    return sstm.str();
}

int Coords::isRight(const Coords* c1, const Coords* c2) const {
    double d = (c2->getY() - c1->getY()) * (x - c2->getX()) - (c2->getX() - c1->getX()) * (y - c2->getY());

    if (std::abs(d) < 0.000001) {
        return 0;
    }
    if (d < 0) {
        return -1;
    }
    if (d > 0) {
        return 1;
    }

}

bool Coords::operator<(const Coords& c) const {
    return this->getX() < c.getX() or(this->getX() == c.getX() and this->getY() < c.getY());
}

bool Coords::operator==(const Coords& c) const {
    double tol = 0.05;
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

