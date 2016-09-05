//============================================================================
// Name        : lcpc.cpp
// Author      : Elias Annila
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "Coords.h"
#include "Funnel.h"
#include "minHeap.h"
#include "lcpfinder.h"
#include "geomfunc.h"
#include "../lib/poly2tri.h"
#include <vector>
#include <gdal/gdal.h>
#include <gdal_priv.h>
#include <gdal/ogrsf_frmts.h>

std::pair<p2t::Point*,p2t::Point*> readCostSurface(char* costSurface, char* targets, char* start, LcpFinder finder) {
    OGRDataSource *csDS;
    csDS = OGRSFDriverRegistrar::Open(costSurface);
    if (csDS == NULL) {
        std::cout << "COST SURFACE NOT FOUND!\n";
        return -1;
    }
    OGRDataSource *targetDS;
    targetDS = OGRSFDriverRegistrar::Open(targets);
    if (targetDS == NULL) {
        std::cout << "TARGET POINTS NOT FOUND!\n";
        return -1;
    }
    OGRDataSource *startDS;
    startDS = OGRSFDriverRegistrar::Open(start);
    if (startDS == NULL) {
        std::cout << "START POINT NOT FOUND!\n";
        return -1;
    }

    OGRLayer* csLr = csDS->GetLayer(0);
    OGRLayer* targetLr = targetDS->GetLayer(0);
    OGRLayer* startLr = startDS->GetLayer(0);
    std::cout << csLr->GetFeatureCount() << " cost surface features found" << std::endl;
    std::cout << targetLr->GetFeatureCount() << " target points found" << std::endl;
    std::cout << startLr->GetFeatureCount() << " start point found" << std::endl;
    csLr->ResetReading();
    targetLr->ResetReading();
    startLr->ResetReading();
    OGRFeature *csFtre;
    OGRFeature *targetFtre;
    OGRFeature *startFtre = startLr->GetNextFeature();
    OGRGeometry* startGeom = startFtre->GetGeometryRef();
    OGRPoint* startPoint;
    p2t::Point* startp2t;
    std::vector<p2t::Point*> targetp2t;
    if (startGeom != NULL
            && wkbFlatten(startGeom->getGeometryType()) == wkbPoint) {
        startPoint = (OGRPoint*) startGeom;
        startp2t = new p2t::Point(startPoint->getX(), startPoint->getY());
    } else {
        printf("no start geometry\n");
    }
    int pIdx = 0;
    while ((csFtre = csLr->GetNextFeature()) != NULL) {
        OGRGeometry* csGeometry = csFtre->GetGeometryRef();
        if (csGeometry != NULL
                && wkbFlatten(csGeometry->getGeometryType()) == wkbPolygon) {
            OGRPolygon *csPolygon = (OGRPolygon *) csGeometry;
            OGRLinearRing* extRing = csPolygon->getExteriorRing();

            std::vector<std::vector < p2t::Point*>> polygon;
            polygon.push_back(std::vector<p2t::Point*>{});

            for (unsigned int i = 0; i < extRing->getNumPoints(); i++) {
                p2t::Point* point = new p2t::Point(extRing->getX(i), extRing->getY(i));
                polygon[0].push_back(point);
            }

            for (unsigned int ri = 0; ri < csPolygon->getNumInteriorRings(); ri++) {
                polygon.push_back(std::vector<p2t::Point*>{});
                OGRLinearRing* intRing = csPolygon->getInteriorRing(ri);
                for (unsigned int i = 0; i < intRing->getNumPoints(); i++) {
                    p2t::Point* point = new p2t::Point(extRing->getX(i), extRing->getY(i));
                    polygon[ri + 1].push_back(point);
                }
                intermidiatePoints(&polygon, 1000);
                finder.addPolygon(polygon, 1);
            }

            if (csPolygon->IsPointOnSurface(startPoint)) {
                finder.addSteinerPoint(startp2t, pIdx);
            }
            
            // ADD TARGET POINTS TO TARGETS!
            pIdx++;
        } else {
            printf("no polygon geometry\n");
        }
    }
    OGRDataSource::DestroyDataSource(csDS);
    return std::pair<p2t::Point*,std::vector<p2t::Point*>(startp2t, targetp2t);
}

int main() {
    /*
     * HOW TO MAKE FAST:
     * A*
     * Data structures: Coordsin naapurimappi, Coordmappi arrayksi? Coords luokka kokonaan pois?
     * Minheap, parempi update operaatio
     *
     *TODO Integrate to QGIS plugin
     *TODO Implement A*
     *TODO Writing results to shapefile (in QGIS)
     *TODO Weakly simple polygon splitting
     *TODO Performance...
     *
     */

    OGRRegisterAll();
    LcpFinder finder{};
    p2t::Point* sp = readCostSurface("testarea.shp","targets.shp","start.shp", finder);



    /*

    std::vector<p2t::Point*> stp = {new p2t::Point
        {0.5, 0.5}};


    finder.addSteinerPoints(stp, 0);
    std::cout << "Steiner added\n";
    std::tr1::unordered_set<Coords, CoordsHasher, std::equal_to<Coords>, std::allocator<Coords> > cmap = finder.getCoordmap();
    for (Coords c : cmap) {
        std::cout << c.toString() << std::endl;
    }
     */
    std::cout << "HALFWAY\n";
    //finder.addPolygon(1, p2, 1, 0.5);

    std::vector<Coords> targets{Coords(1, 2)};
    targets = finder.leastCostPath(Coords(), targets);
    std::cout << "lcp done\n";
    for (Coords goal : targets) {
        while (goal.getPred() != 0) {
            std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;
            goal = *goal.getPred();
        }
        std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;

    }

    return 0;
}

