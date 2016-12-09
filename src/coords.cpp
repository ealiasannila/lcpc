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
    this->target = false;
    this->linePt = false;
}

Coords::Coords(double newx, double newy, int polygon, bool target) {
    x = newx;
    y = newy;
    toStart = -1;
    this->addToPolygon(polygon);
    this->predecessor = 0;
    this->target = target;
    this->linePt = false;
}

Coords::Coords(double newx, double newy, int polygon, bool target, bool linePt) {
    x = newx;
    y = newy;
    toStart = -1;
    this->addToPolygon(polygon);
    this->predecessor = 0;
    this->target = target;
    this->linePt = linePt;
}

Coords::Coords() {
    toStart = -1;
    x = -1;
    y = -1;
    this->predecessor = 0;
    this->target = false;
    this->linePt = false;
}

std::vector<std::pair<const Coords*, double>>*Coords::getNeighbours(int polygon) const {
    try {
        return &(this->neighbours.at(polygon));
    } catch (const std::out_of_range& oor) {
        std::cout << std::fixed << this->getX() << "," << this->getY() << std::endl;
        std::cout << "asking neighbours from polygon where doesn't belong\n";
        std::cout << "getting from: " << polygon << std::endl;
        std::cout << "belongs to:";
        for (int p : this->belongsToPolygons()) {
            std::cout << " " << p;
        }
        std::cout << std::endl;
        exit(1);
    }
}

std::vector<int> Coords::belongsToPolygons() const {
    std::vector<int> polygons;
    for (std::pair<int, std::vector<std::pair<const Coords*, double>>> p : this->neighbours) {
        polygons.push_back(p.first);
    }
    return polygons;
}

bool Coords::belongsToPolygon(int p) const {
    return this->neighbours.find(p) != this->neighbours.end();
}

void Coords::addToPolygon(int polygon) const {

    this->neighbours.insert(std::make_pair(polygon, std::vector<std::pair<const Coords*, double>>
    {
    }));

}

void Coords::addNeighbours(const Coords* c, int polygon, double friction) const {
    if (c->linePt) {
        
        std::cout<<std::fixed<<c->getX()<<","<<c->getY()<<std::endl;
        std::cout << friction << std::endl;
    }
    try {

        auto it = std::lower_bound(neighbours.at(polygon).begin(), neighbours.at(polygon).end(), c, std::bind(&Coords::compareNeighbours, this, std::placeholders::_1, std::placeholders::_2, polygon));
        if ((it == neighbours.at(polygon).end() or it->first != c) and (neighbours.at(polygon).empty() or c != neighbours.at(polygon).at(0).first)) {
            neighbours.at(polygon).insert(it, std::make_pair(c, friction));
        } else if (it->first == c and it->second > friction) {
            std::cout << "updating friction\n";
            it->second = friction;
        }
    } catch (const std::out_of_range& oor) {
        std::cout << "trying to add neighbours to polygon where doesn't belong\n";
        exit(1);
    }
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

