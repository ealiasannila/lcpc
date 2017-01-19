/*
 * lcpfinder.h
 *
 *  Created on: Aug 25, 2016
 *      Author: elias
 */

#ifndef SRC_LCPFINDER_H_
#define SRC_LCPFINDER_H_

#include "coords.h"
#include "funnel.h"
#include "min_heap.h"
#include "geomfunc.h"
#include "../lib/poly2tri.h"
#include <tr1/unordered_set>
#include "defs.h"
#include<limits>
#include "compare.h"
#include <forward_list>

struct CoordsHasher {

    std::size_t operator()(const Coords& c) const {
        std::size_t res = 17;
        res = res * 31 + std::tr1::hash<double>()((int) c.getX());
        res = res * 31 + std::tr1::hash<double>()((int) c.getY());
        return res;
    }
};

class LcpFinder {
private:
    int numOfTargets = 0;
    double maxDist = 0.0;
    double minFriction = std::numeric_limits<double>::infinity();
    std::tr1::unordered_set<Coords, CoordsHasher> coordmap;
    std::vector<std::vector<std::vector<const Coords*>>> polygons;
    std::map<int, std::forward_list<const Coords*>> targetPoints;
    std::map<int, std::forward_list<const Coords*>> linePoints;
    std::map<int, std::vector<const Coords*>> fences; //first point of a fence in polygon (either at border or first point of fence)
    std::vector<double> frictions;
    std::vector<bool> triangulated;

    std::deque<Funnel> initFQue(const Coords* c, int polygon, nSet*neighbours);
    void findNeighboursInPolygon(const Coords* c, int polygon, nSet* neighbours);
    double toClosestEnd(const Coords* c);
    void checkTargets(int polygon, Triangle* newTri);
    void checkLinear(int polygon, Triangle* newTri, bool subtri);
    double checkFences(int polygon,  const Coords* a, const Coords* b);
    void subTriangles(Triangle* newTri, int polygon, const Coords* c);    
    std::vector<const Coords*> segmentPolygonIntersection(int polygon,const Coords* a, const Coords* b, int* next);
    
public:
    void addBuffers( double d);
    void addBufferPoint(const Coords* prev, const Coords* c, const Coords* next, double d,  int polygon);
    void setMaxD(double d);
    nSet findNeighbours(const Coords* c);

    int getPolygonCount() {
        return this->polygons.size();
    };

    std::tr1::unordered_set<Coords, CoordsHasher>* getCoordmap() {
        return &coordmap;
    }
    std::deque<const Coords*> leastCostPath(int algorithm);
    void addPolygon(std::vector<std::vector<p2t::Point*>> points, double friction);
    //void addLine(std::vector<p2t::Point*>*, double frictionForwards,double frictionBackwards, std::array<int,2> polygons);
    void addStartPoint(p2t::Point* start, int polygon);
    void addTargetPoint(p2t::Point* steinerpoint, int polygon);
    const Coords* addLinePoint(p2t::Point* linepoint, int polygon);
    const Coords* addLinePoint(std::array<double,2>, int polygon);
    void addLine(std::vector<p2t::Point*>* points, double friction, double crossing);
    void triangulate(int polygon);
    const Coords* startPoint;
    std::array<int, 2> containingPolygon(p2t::Point* p);

    ~LcpFinder();
};



#endif /* SRC_LCPFINDER_H_ */
