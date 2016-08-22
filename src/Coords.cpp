/*
 * Coords.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#include "Coords.h"
#include <sstream>
#include <math.h>

Coords::Coords(double newx, double newy) {
	x = newx;
	y = newy;
}

Coords::Coords(double newx, double newy, int polygon) {
	x = newx;
	y = newy;
	this->leftNeighbours.insert(std::pair<int, std::set<Coords*>>(polygon, std::set<Coords*>()));
	this->rightNeighbours.insert(std::pair<int, std::set<Coords*>>(polygon, std::set<Coords*>()));
}

Coords::Coords(){
	x = -1;
	y = -1;
}

std::pair<std::set<Coords*>::iterator,std::set<Coords*>::iterator> Coords::getRightNeighbours(int polygon){
	return std::pair<std::set<Coords*>::iterator,std::set<Coords*>::iterator>(this->rightNeighbours.at(polygon).begin(),this->rightNeighbours.at(polygon).end())  ;
}
std::pair<std::set<Coords*>::iterator,std::set<Coords*>::iterator> Coords::getLeftNeighbours(int polygon){
	return std::pair<std::set<Coords*>::iterator,std::set<Coords*>::iterator>(this->leftNeighbours.at(polygon).begin(),this->leftNeighbours.at(polygon).end())  ;
}

std::pair<std::map<int,std::set<Coords*>>::iterator,std::map<int,std::set<Coords*>>::iterator> Coords::belongsToPolygons(){
	return std::pair<std::map<int,std::set<Coords*>>::iterator,std::map<int,std::set<Coords*>>::iterator>(this->leftNeighbours.begin(),this->leftNeighbours.end());
}

void Coords::addToPolygon(int polygon){
	this->leftNeighbours.insert(std::pair<int, std::set<Coords*>>(polygon, std::set<Coords*>()));
	this->rightNeighbours.insert(std::pair<int, std::set<Coords*>>(polygon, std::set<Coords*>()));
}

void Coords::addNeighbours(Coords* l, Coords* r, int polygon){
	leftNeighbours.at(polygon).insert(l);
	rightNeighbours.at(polygon).insert(r);
}

std::string Coords::toString(){
	std::stringstream sstm;
	sstm <<"X: "<<x<<" Y: "<<y;
	return  sstm.str();
}
int Coords::isRight(Coords* c1, Coords* c2){
	double d = (c2->getY() - c1->getY()) * (x - c2->getX())
	                - (c2->getX() - c1->getX()) * (y - c2->getY());
	 if(d<0){
	 		 return -1;
	 }if(d>0){
		 return 1;
	 }
	 return 0;
}

double Coords::eucDistSquared(Coords* c1) {
	return pow(c1->getX()-x, 2) + pow(c1->getY()-y,2);
}

double Coords::eucDist(Coords* c1){
	return sqrt(this->eucDistSquared(c1));
}

Coords::~Coords() {
	// TODO Auto-generated destructor stub
}

