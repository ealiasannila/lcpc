/*
 * Funnel.h
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */
#ifndef FUNNEL_H_
#define FUNNEL_H_

#include <deque>
#include "coords.h"
#include <string>


class Funnel {
public:
	const Coords* firstLeft;
	const Coords* firstRight;
	const Coords* apex;
        const Triangle* t;
        int base;
        int inFirstSector(const Coords* o);
	int getNewBase(const Triangle* nt, const Coords* oldOpposing, int LeftOrRight);

	Funnel(const Coords* a, const Triangle* t);
	Funnel(const Coords* fl, const Coords* fr, const Coords* a, const Triangle* t, int b);
        const Coords* getOpposing();
	Funnel split();
	void stepForward(std::deque<Funnel>* funnelQueue, nSet* neighbours, double friction);
	std::string toString();
	const Coords* getApex(){return this->apex;}
	virtual ~Funnel();
};

#endif /* FUNNEL_H_ */
