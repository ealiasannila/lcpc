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

Coords::Coords(double newx, double newy, int polygon, int id) {
    x = newx;
    y = newy;
    toStart = -1;
    this->addToPolygon(polygon);
    this->predecessor = 0;
    this->id = id;
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

nContainer* Coords::getNeighbours(int polygon) const {
    try {
        return &this->neighbours.at(polygon);
    } catch (const std::out_of_range& oor) {
        std::cout << "NoNeighbours: \n";
        std::cout << "polygon: " << polygon;
        exit(1);
    }
}

std::vector<int> Coords::belongsToPolygons() const {
    std::vector<int> polygons;
    for (std::pair<int, nContainer> p : this->leftNeighbours) {
        polygons.push_back(p.first);
    }
    return polygons;
}

void Coords::addToPolygon(int polygon) const {
    this->leftNeighbours.insert(std::pair<int, nContainer>(polygon, nContainer()));
    this->rightNeighbours.insert(std::pair<int, nContainer>(polygon, nContainer()));
    this->neighbours.insert(std::pair<int, nContainer>(polygon, nContainer()));
}

void Coords::addNeighbours(const Coords* l, const Coords* r, int polygon) const {
    try {

        auto li = std::lower_bound(neighbours.at(polygon).begin(), neighbours.at(polygon).end(), l);

        if (neighbours.at(polygon).empty() or *(li) != l) {
            neighbours.at(polygon).insert(li, l);
        }
        
        auto ri = std::lower_bound(neighbours.at(polygon).begin(), neighbours.at(polygon).end(), r);

        if (neighbours.at(polygon).empty() or(*ri) != r) {
            neighbours.at(polygon).insert(ri, r);
        }


        leftNeighbours.at(polygon).push_back(l);
        rightNeighbours.at(polygon).push_back(r);
    } catch (const std::out_of_range& oor) {
        std::cout << "NoNeighbours: \n";
        std::cout << this->toString();
        std::cout << "polygon: " << polygon;
        exit(1);
    }
}

void Coords::setToStart(double cost) const {
    toStart = cost;
}

void Coords::setPred(const Coords* pred) const {
    predecessor = pred;
}

std::string Coords::toString() const {
    std::stringstream sstm;
    sstm << "ID: " << id; // << "X: " << x << " Y: " << y << " Polygons: ";
    /*
    for (std::pair<int, nContainer> p : this->leftNeighbours) {
        sstm << p.first << " ";
    }*/
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

