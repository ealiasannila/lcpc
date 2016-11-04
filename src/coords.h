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
#include "sorted_vector.h"

class Coords {
private:

    mutable double toStart;
    mutable double toEnd;
    mutable const Coords* predecessor;
    mutable std::tr1::unordered_map<int, SortedVector<const Coords*>> neighbours;
    mutable allNContainer leftNeighbours;
    mutable allNContainer rightNeighbours;
    double x;
    double y;


public:
    bool target;
    Coords();
    Coords(double x, double y);
    Coords(double x, double y, int polygon, bool target);

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
    int isRight(const Coords* c1, const Coords* c2) const;
    nContainer* getRightNeighbours(int polygon) const;
    nContainer* getLeftNeighbours(int polygon) const;
    SortedVector<const Coords*>* getNeighbours(int polygon) const;
    std::vector<int> belongsToPolygons() const;
    void addToPolygon(int polygon) const;
    void addNeighbours(const Coords* l, const Coords* r, int polygon) const;
    std::string toString() const;
    virtual ~Coords();


    bool operator<(const Coords& c) const;

    bool operator==(const Coords& c) const;
    bool operator!=(const Coords& c) const;
};

#endif /* const Coords_H_ */
