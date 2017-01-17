/*
 * const Coords.h
 *
 *  Created on: Aug 16, 2016
 *      Author: elias
 */

#ifndef  Coords_H_
#define  Coords_H_
#include <string>
#include "defs.h"
#include <set>
#include<iostream>
#include<boost/heap/fibonacci_heap.hpp>

#include "compare.h"
class Coords;

struct Triangle {
    std::array<const Coords*, 3> points;
    std::array<const Triangle*, 3> neighbours;
    std::vector<const Coords*> interiorPoints;
};

class Coords {
private:
    mutable double toStart;
    mutable double toEnd;
    mutable const Coords* predecessor;
    mutable std::map<int, std::vector<const Triangle*>> triangles;
    double x;
    double y;


public:
    /*
     * -1 temporary point
     *  0 normal point
     *  1 target point
     *  2 linear point
     */
    int flag = 0;
    mutable std::vector<std::tuple<const Coords*, double, double>> linearNeighbours;

    struct cmpr {

        /*
        bool (*f)(const Coords*, const Coords*);

        cmpr() {
            f = &compAstar;
        }

        cmpr(bool (*comparefunction)(const Coords*, const Coords*)) {
            f = comparefunction;
        }
         */
        bool operator()(const Coords* x, const Coords* y) const {
            //return (x->getToStart() + x->getToEnd() > y->getToStart() + y->getToEnd());
            return (x->getToStart() > y->getToStart());
        }
    };
    mutable boost::heap::fibonacci_heap<const Coords*, boost::heap::compare<cmpr>>::handle_type handle;
    Coords();
    Coords(double x, double y);
    Coords(double x, double y, int polygon);
    Coords(double x, double y, int polygon, int flag);

    double getX() const {
        return x;

    }

    double getToStart() const {
        return toStart;
    }

    double getToEnd() const {
        return toEnd;
    }

    const Coords* getPred() const {
        return predecessor;
    }
    void setToStart(double cost) const;
    void setToEnd(double cost) const;
    void setPred(const Coords* pred) const;

    double getY() const {
        return y;
    }
    bool belongsToPolygon(int p) const;
    int isRight(const Coords* c1, const Coords* c2) const;
    std::vector<const Triangle*>*getTriangles(int polygon) const;
    void addTriangle(Triangle* t, int p) const;
    std::vector<int> belongsToPolygons() const;
    void addToPolygon(int polygon) const;
    void addLinearNeighbour(const Coords* n, double friction) const;
    void addLinearNeighbour(const Coords* n, double friction, double crossingcost) const;
    std::string toString() const;
    virtual ~Coords();


    bool operator<(const Coords& c) const;

    bool operator==(const Coords& c) const;
    bool operator!=(const Coords& c) const;
};



#endif /* const Coords_H_ */
