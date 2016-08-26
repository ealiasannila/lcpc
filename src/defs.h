/*
 * defs.h
 *
 *  Created on: Aug 24, 2016
 *      Author: elias
 */

#ifndef SRC_DEFS_H_
#define SRC_DEFS_H_

#include <map>
#include <list>
#include <queue>
#include <set>
#include <vector>
#include <array>

class Coords;
typedef std::vector<const Coords*> nContainer;
typedef std::map<int, nContainer> allNContainer;
typedef std::pair<nContainer::iterator, nContainer::iterator>	neighIter;
typedef std::pair<std::map<int, nContainer>::iterator, std::map<int, nContainer>::iterator> allNeighIter;
typedef std::map<const Coords*, int> nSet; //used in dijsktra/A* must not allow duplicates.
#endif /* SRC_DEFS_H_ */
