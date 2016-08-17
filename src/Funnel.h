/*
 * Funnel.h
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */
#ifndef FUNNEL_H_
#define FUNNEL_H_

#include <deque>
#include "Coords.h"
#include "Edge.h"

class Funnel {
private:
	Coords* apex;
	std::deque<Coords*> sides;

public:
	Funnel(Coords* l, Coords* a, Coords* r);
	Edge getBase();
	void addR(Coords* c);
	void addL(Coords* c);
	Coords* getApex(){return apex;}
	virtual ~Funnel();
};

#endif /* FUNNEL_H_ */
