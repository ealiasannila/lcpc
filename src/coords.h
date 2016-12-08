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
struct RotaryCompare;
class Coords {
private:

    mutable double toStart;
    mutable double toEnd;
    mutable const Coords* predecessor;
    mutable std::map<int, std::vector<std::pair<const Coords*, double>>> neighbours;
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
    std::vector<std::pair<const Coords*, double>>*getNeighbours(int polygon) const;
    std::vector<int> belongsToPolygons() const;
    void addToPolygon(int polygon) const;
    void addNeighbours(const Coords* c, int polygon, double friction) const;
    std::string toString() const;
    virtual ~Coords();


    bool operator<(const Coords& c) const;

    bool operator==(const Coords& c) const;
    bool operator!=(const Coords& c) const;
};

struct RoataryCompare {
    const Coords* center;

    RoataryCompare(const Coords* c) {
        center = c;
    }
    bool operator()( std::pair<const Coords*, double> p,const Coords* c) const {
        return center->isRight(p.first, c)==1;
    }
};

#endif /* const Coords_H_ */
