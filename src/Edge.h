/*
 * Edge.h
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#ifndef EDGE_H_
#define EDGE_H_

#include "Coords.h"

class Edge {
private:
	Coords* l;
	Coords* r;
public:
	Edge(Coords* l, Coords* r);
	Coords* getL(){return l;}
	Coords* getR(){return r;}
	virtual ~Edge();
};

#endif /* EDGE_H_ */
