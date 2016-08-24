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

nContainer Coords::getRightNeighbours(int polygon) const {
	return this->rightNeighbours.at(polygon);
}
nContainer Coords::getLeftNeighbours(int polygon) const {
	return this->leftNeighbours.at(polygon);
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
	leftNeighbours.at(polygon).push_back(l);
	rightNeighbours.at(polygon).push_back(r);
}

void Coords::setToStart(double cost) const{
	toStart = cost;
}
void Coords::setPred(const Coords* pred) const{
	predecessor = pred;
}

std::string Coords::toString() const {
	std::stringstream sstm;
	sstm << "X: " << x << " Y: " << y;
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

double Coords::eucDistSquared(const Coords* c1) const{
	return pow(c1->getX() - x, 2) + pow(c1->getY() - y, 2);
}

double Coords::eucDist(const Coords* c1) const {
	return sqrt(this->eucDistSquared(c1));
}

bool Coords::operator==(const Coords& c) const {
	double tol = 0.000001;
	if (std::abs(c.getX() - this->getX()) < tol and std::abs(c.getY() - this->getY()) < tol) {
		return true;
	}
	return false;
}

bool Coords::operator!=(const Coords& c) const {
	return !this->operator ==(c);
}

Coords::~Coords() {
	// TODO Auto-generated destructor stub
}

