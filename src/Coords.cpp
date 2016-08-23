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
 Coords:: Coords(double newx, double newy) {
	x = newx;
	y = newy;
}

 Coords:: Coords(double newx, double newy, int polygon) {
	x = newx;
	y = newy;
	this->leftNeighbours.insert(std::pair<int, std::set< const Coords*>>(polygon, std::set< const Coords*>()));
	this->rightNeighbours.insert(std::pair<int, std::set< const Coords*>>(polygon, std::set< const Coords*>()));
}

 Coords:: Coords() {
	x = -1;
	y = -1;
}

std::pair<std::set< const Coords*>::iterator, std::set< const Coords*>::iterator>  Coords::getRightNeighbours(int polygon) const {
	return std::pair<std::set< const Coords*>::iterator, std::set< const Coords*>::iterator>(this->rightNeighbours.at(polygon).begin(), this->rightNeighbours.at(polygon).end());
}
std::pair<std::set< const Coords*>::iterator, std::set< const Coords*>::iterator>  Coords::getLeftNeighbours(int polygon)const {
	return std::pair<std::set< const Coords*>::iterator, std::set< const Coords*>::iterator>(this->leftNeighbours.at(polygon).begin(), this->leftNeighbours.at(polygon).end());
}

std::pair<std::map<int, std::set< const Coords*>>::iterator, std::map<int, std::set< const Coords*>>::iterator>  Coords::getAllLeftN() const{
	return std::pair<std::map<int, std::set< const Coords*>>::iterator, std::map<int, std::set< const Coords*>>::iterator>(this->leftNeighbours.begin(), this->leftNeighbours.end()) ;
}
std::pair<std::map<int, std::set< const Coords*>>::iterator, std::map<int, std::set< const Coords*>>::iterator>  Coords::getAllRightN() const{
	return std::pair<std::map<int, std::set< const Coords*>>::iterator, std::map<int, std::set< const Coords*>>::iterator>(this->rightNeighbours.begin(), this->rightNeighbours.end());
}

void Coords::addToPolygon(int polygon)const {
	this->leftNeighbours.insert(std::pair<int, std::set< const Coords*>>(polygon, std::set< const Coords*>()));
	this->rightNeighbours.insert(std::pair<int, std::set< const Coords*>>(polygon, std::set< const Coords*>()));
}

void Coords::addNeighbours(const Coords* l, const Coords* r, int polygon) const {
	leftNeighbours.at(polygon).insert(l);
	rightNeighbours.at(polygon).insert(r);
}

std::string  Coords::toString() const {
	std::stringstream sstm;
	sstm << "X: " << x << " Y: " << y;
	return sstm.str();
}
int  Coords::isRight(const Coords* c1, const Coords* c2) const {
	double d = (c2->getY() - c1->getY()) * (x - c2->getX()) - (c2->getX() - c1->getX()) * (y - c2->getY());
	if (d < 0) {
		return -1;
	}
	if (d > 0) {
		return 1;
	}
	return 0;
}

double  Coords::eucDistSquared(const Coords* c1) {
	return pow(c1->getX() - x, 2) + pow(c1->getY() - y, 2);
}

double  Coords::eucDist(const Coords* c1) {
	return sqrt(this->eucDistSquared(c1));
}


bool Coords::operator==( const Coords& c)  const{
	double tol = 0.000001;
	if (std::abs(c.getX() - this->getX()) < tol and std::abs(c.getY() - this->getY()) < tol) {
		return true;
	}
	return false;
}

bool Coords::operator!=( const Coords& c) const {
	return !this->operator ==(c);
}

 Coords::~ Coords() {
	// TODO Auto-generated destructor stub
}

