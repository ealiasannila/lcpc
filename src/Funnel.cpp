/*
 * Funnel.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#include "Funnel.h"

Funnel::Funnel(Coords* l, Coords* a, Coords* r) {
	apex = a;
	sides.push_back(r);
	sides.push_front(l);
}

//TODO varaudu jos deque on tyhjä tai jos törmää apexiin
Edge Funnel::getBase(){
	return Edge(sides.front(),sides.back());
}

void Funnel::addR(Coords* c){
	sides.push_back(c);
}
void Funnel::addL(Coords* c){
	sides.push_front(c);
}

Funnel::~Funnel() {
	// TODO Auto-generated destructor stub
}


