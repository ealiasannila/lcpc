/*
 * helperfunctions.hpp
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */


#ifndef HELPERFUNCTIONS_OMA_HPP
#define HELPERFUNCTIONS_OMA_HPP

#include "Coords.h"
#include <math.h>

double eucDistSquared(Coords c1, Coords c2) {
	return pow(c1.getX()-c2.getX(), 2) + pow(c1.getY()-c2.getY(),2);
}

double eucDist(Coords c1, Coords c2){
	return sqrt(eucDistSquared(c1,c2));
}

/*returns -1 if c3 is left of line between c1 and c2,
1 if right of the line and 0 if all points are on same line*/
int isRight(Coords c1, Coords c2, Coords c3){
	 double d = (c2.getY() - c1.getY()) * (c3.getX() - c2.getX())
	                - (c2.getX() - c1.getX()) * (c3.getY() - c2.getY());
	 if(d<0){
	 		 return -1;
	 }if(d>0){
		 return 1;
	 }
	 return 0;
}

#endif /* HELPERFUNCTIONS_HPP_ */
