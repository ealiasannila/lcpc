/*
 * defs.h
 *
 *  Created on: Aug 24, 2016
 *      Author: elias
 */

#ifndef SRC_DEFS_H_
#define SRC_DEFS_H_

#include <tr1/unordered_map>
#include <map>
#include <vector>


class Coords;
typedef std::vector<const Coords*> nContainer;
typedef std::tr1::unordered_map<int, nContainer> allNContainer;
typedef std::map<const Coords*, double> nSet; //used in dijsktra/A* must not allow duplicates.
#endif /* SRC_DEFS_H_ */
