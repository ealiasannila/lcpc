//============================================================================
// Name        : lcpc.cpp
// Author      : Elias Annila
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "tests/test.h"
#include "../lib/geometry.h"
#include <algorithm>

Coords* getOpposing(Coords* l, Coords* r){
	std::vector<Coords*> cn(2); //common neighbours (may be 1 or 2)
	std::vector<Coords*>::iterator it = std::set_intersection(l->getLeftNeighbours(), l->getLeftNeighboursEnd(), r->getRightNeighbours(), r->getRightNeighboursEnd(), cn.begin());
	cn.resize(it-cn.begin());
	for (it=cn.begin(); it!=cn.end(); ++it){
		Coords* c = *it;
		if(c->isRight(l, r)==-1){
			return c;
		}
	}

	std::cout<<"No opposing found"<<std::endl;
	return 0;
}

int main() {
	//test();
	Polygon poly("sample1.bdm", true);
	poly.setDebugOption(false);      //set debug flag;
	poly.triangulation();
	PointbaseMap points = poly.points();
	Triangles triangles = poly.triangles();



	for (std::list<Triangle>::iterator it = triangles.begin();it != triangles.end(); it++) {
		Triangle triangle = *it;
		std::cout<<std::endl;
		Coords c [3];
		for(unsigned i = 0; i<3; i++){
			Pointbase pb = *points.at(triangle[i]);
			c[i] = Coords(pb.x, pb.y);
		}
		// if triangle orientation is clockwise turn it to CCW
		if(c[0].isRight(&c[1], &c[2])==1){
			Coords tmp = c[2];
			c[2] = c[1];
			c[1] = tmp;
		}
		c[0].addNeighbours(&c[2], &c[1]);
		c[1].addNeighbours(&c[0], &c[2]);
		c[2].addNeighbours(&c[1], &c[0]);
		std::cout<<"--"<<std::endl;
		Coords* o = getOpposing(&c[0],&c[2]);
		if(o != 0){
			std::cout<<o->getX()<< "-" <<o->getY()<<std::endl;
		}
	}
	return 0;
}


