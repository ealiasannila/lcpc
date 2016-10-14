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
#include <string>


class Funnel {
private:
	std::vector<const Coords*> lc;
	std::vector<const Coords*> rc;
	void shrink(const Coords* o, std::vector<const Coords*>* chain, unsigned lastRemaining);
	int inFirstSector(const Coords* o);
	int findInChain(const Coords* o, std::vector<const Coords*> chain, int side);

public:
	Funnel(const Coords* l, const Coords* a, const Coords* r);
	Funnel(std::vector<const Coords*> lc, std::vector<const Coords*> rc);
	std::pair<const Coords*, const Coords*> getBase();
	Funnel split(const Coords* o);
	void reactToOpposite(const Coords* o, std::deque<Funnel>* funnelQueue, nSet* neighbours, int polygon);
	std::vector<const Coords*> getLC();
	std::vector<const Coords*> getRC();
	std::string toString();
	const Coords* getApex(){return lc.front();}
	virtual ~Funnel();
};

#endif /* FUNNEL_H_ */
