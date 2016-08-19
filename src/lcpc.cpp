//============================================================================
// Name        : lcpc.cpp
// Author      : Elias Annila
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "tests/test.h"
#include "../lib/geometry.h"



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
		if(c[0].isRight(&c[1], &c[2])==1){
			Coords tmp = c[2];
			c[2] = c[1];
			c[1] = tmp;
		}
		c[0].addNeighbours(&c[2], &c[1]);
		c[1].addNeighbours(&c[0], &c[2]);
		c[2].addNeighbours(&c[1], &c[0]);


	}
	return 0;
}


