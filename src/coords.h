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
    bool linePt;
    Coords();
    Coords(double x, double y);
    Coords(double x, double y, int polygon, bool target);
    Coords(double x, double y, int polygon, bool target, bool linePt);

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
    std::vector<std::pair<const Coords*, double>>*getNeighbours(int polygon) const;
    std::vector<int> belongsToPolygons() const;
    void addToPolygon(int polygon) const;
    void addNeighbours(const Coords* c, int polygon, double friction) const;
    std::string toString() const;
    virtual ~Coords();

    bool compareNeighbours(std::pair<const Coords*, double> p, const Coords* c, int polygon) const {
        const Coords* first;
        try {
            first = this->neighbours.at(polygon).at(0).first;
        } catch (const std::out_of_range& oor) {
            return this->isRight(p.first, c) == 1;
        }
        if(p.first == first){
            return true;
        }if(c == first){
            return true;
        }
        int p_half = p.first->isRight(first, this);
        int c_half = c->isRight(first, this);
       
        if(p_half != c_half){
            return p_half<c_half;
        }
        return this->isRight(p.first, c) == 1;


    }
    bool operator<(const Coords& c) const;

    bool operator==(const Coords& c) const;
    bool operator!=(const Coords& c) const;
};



#endif /* const Coords_H_ */
