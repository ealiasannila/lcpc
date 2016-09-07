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

void readCostSurface(const char* costSurface, const char* targets, const char* start, LcpFinder* finder) {
    OGRDataSource *csDS;
    csDS = OGRSFDriverRegistrar::Open(costSurface);
    if (csDS == NULL) {
        std::cout << "COST SURFACE NOT FOUND!\n";
        return;
    }
    OGRDataSource *targetDS;
    targetDS = OGRSFDriverRegistrar::Open(targets);
    if (targetDS == NULL) {
        std::cout << "TARGET POINTS NOT FOUND!\n";
        return;
    }
    OGRDataSource *startDS;
    startDS = OGRSFDriverRegistrar::Open(start);
    if (startDS == NULL) {
        std::cout << "START POINT NOT FOUND!\n";
        return;
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
    std::vector<OGRPoint*> targetOGRPoints;

    while ((targetFtre = targetLr->GetNextFeature()) != NULL) {
        OGRGeometry* targetGeom = targetFtre->GetGeometryRef();
        if (targetGeom != NULL
                && wkbFlatten(targetGeom->getGeometryType()) == wkbPoint) {
            OGRPoint* targetPoint = (OGRPoint*) targetGeom;
            targetOGRPoints.push_back(targetPoint);
        } else {
            printf("no target geometry\n");
        }

    }
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
            if (extRing->isClockwise()) {
                extRing->reverseWindingOrder();
            }

            std::vector<std::vector < p2t::Point*>> polygon;
            polygon.push_back(std::vector<p2t::Point*>{});

            int ringsize = extRing->getNumPoints();


            for (unsigned int i = 0; i < ringsize; i++) {
                p2t::Point* point = new p2t::Point(extRing->getX(i), extRing->getY(i));
                polygon[0].push_back(point);
                //std::cout<<"adding: "<<point->x<<","<<point->y<<std::endl;
            }
            polygon[0].pop_back();

            for (unsigned int ri = 0; ri < csPolygon->getNumInteriorRings(); ri++) {
                polygon.push_back(std::vector<p2t::Point*>{});
                OGRLinearRing* intRing = csPolygon->getInteriorRing(ri);
                if (intRing->isClockwise()) {
                    intRing->reverseWindingOrder();
                }
                ringsize = intRing->getNumPoints();
                for (unsigned int i = 0; i < ringsize; i++) {
                    p2t::Point* point = new p2t::Point(extRing->getX(i), extRing->getY(i));
                    polygon[ri + 1].push_back(point);
                }
                polygon[ri + 1].pop_back();
            }


            intermidiatePoints(&polygon, 1000);
            finder->addPolygon(polygon, 1);

            if (csPolygon->IsPointOnSurface(startPoint)) {
                finder->addStartPoint(startp2t, pIdx);
                std::cout << "ADDING START POINT\n";
            }
            for (int i = 0; i < targetOGRPoints.size(); i++) {
                if (csPolygon->IsPointOnSurface(targetOGRPoints[i])) {
                    finder->addSteinerPoint(new p2t::Point(targetOGRPoints[i]->getX(), targetOGRPoints[i]->getY()), pIdx);
                }
            }


            // ADD TARGET POINTS TO TARGETS!
            pIdx++;

        } else {
            printf("no polygon geometry\n");
        }
    }

    OGRDataSource::DestroyDataSource(csDS);
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
    readCostSurface("testpolygon.shp", "targets.shp", "start.shp", &finder);
    finder.leastCostPath();

    std::cout << "LCP SEARCH DONE\n";


    std::map<const Coords*, const Coords*> ancestors;
    std::map<const Coords*, std::vector<double>> x;
    std::map<const Coords*, std::vector<double>> y;


    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer("point_out", NULL, wkbLineString, NULL);
    if (poLayer == NULL) {
        printf("Layer creation failed.\n");
        exit(1);
    }

    SHPHandle hSHP = SHPCreate(outFile.toStdString().c_str(), SHPT_ARC);
    std::cout << "handle created\n";


    std::cout << "Initial insert:\n";
    for (const Coords* point : results) {
        ancestors[point] = point;
        x[point] = std::vector<double>{point->getX()};
        y[point] = std::vector<double>{point->getY()};
        std::cout << point->toString() << std::endl;
    }

    std::cout << "Starting loop:\n";
    while (!results.empty()) {
        const Coords* point = results.front();
        const Coords* pred = point->getPred();

        std::cout << "Po: " << point->toString() << std::endl;
        std::cout << "Pr: " << point->toString() << std::endl;


        results.pop_front();
        if (pred == 0) {
            continue;
        }

        const Coords* ancestor = ancestors[point];
        std::cout << "An: " << point->toString() << std::endl;

        if (ancestors.find(pred) != ancestors.end()) {
            ancestors[pred] = pred;
            x[pred] = std::vector<double>{pred->getX()};
            y[pred] = std::vector<double>{pred->getY()};
            std::cout << "Adding point to pred: " << pred->toString() << std::endl;
        } else {
            ancestors[pred] = ancestor;
            results.push_back(pred);
        }
        x[ancestor].push_back(pred->getX());
        y[ancestor].push_back(pred->getY());
        std::cout << "Adding point to ancestor: " << ancestor->toString() << std::endl;

    }
    std::cout << "Paths mapped\n";
    for (std::pair<const Coords*, std::vector<double>> xs : x) {
        std::vector<double> ys = y[xs.first];

        std::cout << "X pointer\n";
        double *xp = &(xs.second[0]);
        std::cout << "Y pointer\n";
        double *yp = &(ys[0]);

        std::cout << "Create SHP\n";
        SHPObject* psObject = SHPCreateSimpleObject(SHPT_ARC, xs.second.size(), xp, yp, NULL);
        std::cout << "Write SHP\n";
        SHPWriteObject(hSHP, -1, psObject);
        SHPDestroyObject(psObject);
    }
    SHPClose(hSHP);

    return 0;
    /*    std::vector<Coords> targets = finder.leastCostPath();
        std::cout << "lcp done\n";
        for (Coords goal : targets) {
            while (goal.getPred() != 0) {
                std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;
                goal = *goal.getPred();
            }
            std::cout << "x: " << goal.getX() << " y: " << goal.getY() << "cost: " << goal.getToStart() << std::endl;

        }
     */
    return 0;
}

