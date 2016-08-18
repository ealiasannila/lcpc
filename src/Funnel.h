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
#include <string>
#include <vector>


class Funnel {
private:
	std::deque<Coords*> lc;
	std::deque<Coords*> rc;
	void shrink(Coords* o, std::deque<Coords*>* chain, unsigned lastRemaining);
	int inFirstSector(Coords* o);
	int findInChain(Coords* o, std::deque<Coords*> chain, int side);

public:
	Funnel(Coords* l, Coords* a, Coords* r);
	Funnel(std::deque<Coords*> lc, std::deque<Coords*> rc);
	Edge getBase();
	Funnel split(Coords* o);
	void reactToOpposite(Coords* o);
	std::deque<Coords*> getLC();
	std::deque<Coords*> getRC();
	std::string toString();
	Coords* getApex(){return lc.front();}
	virtual ~Funnel();
};

#endif /* FUNNEL_H_ */
