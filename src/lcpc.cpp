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


    OGRFeature *startFtre = startLr->GetNextFeature();

    std::vector<OGRPoint*> targetOGRPoints;

    std::vector<OGRFeature*> targetPointers;
    OGRFeature * targetFtre;
    while ((targetFtre = targetLr->GetNextFeature()) != NULL) {
        targetPointers.push_back(targetFtre);
        OGRGeometry* targetGeom = targetFtre->GetGeometryRef();
        if (targetGeom != NULL
                && wkbFlatten(targetGeom->getGeometryType()) == wkbPoint) {
            OGRPoint* targetPoint = (OGRPoint*) targetGeom;
            targetOGRPoints.push_back(targetPoint);
        } else {
            printf("no target geometry\n");
        }

    }


    OGRPoint* startPoint;
    p2t::Point* startp2t;
    OGRGeometry* startGeom = startFtre->GetGeometryRef();

    if (startGeom != NULL
            && wkbFlatten(startGeom->getGeometryType()) == wkbPoint) {
        startPoint = (OGRPoint*) startGeom;
        startp2t = new p2t::Point(startPoint->getX(), startPoint->getY());
    } else {
        printf("no start geometry\n");
    }

    int pIdx = 0;
    OGRFeature * csFtre;
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
            delete polygon[0].back();
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
                delete polygon[ri+1].back();
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
        OGRFeature::DestroyFeature(csFtre);

    }

    for (OGRFeature* targetPt : targetPointers) {
        OGRFeature::DestroyFeature(targetPt);

    }

    OGRFeature::DestroyFeature(startFtre);



    OGRDataSource::DestroyDataSource(csDS);
    OGRDataSource::DestroyDataSource(startDS);
    OGRDataSource::DestroyDataSource(targetDS);


}

void readCostSurfaceDummy(const char* costSurface, const char* targets, const char* start, LcpFinder* finder) {
    std::vector<std::vector < p2t::Point*>> polygon
    {
        {
            new p2t::Point(0, 0), new p2t::Point(1, 0), new p2t::Point(1, 1), new p2t::Point(0, 1)
        }
    };
    finder->addPolygon(polygon, 1);
    finder->addStartPoint(new p2t::Point(0.5, 0.5), 0);
    finder->addSteinerPoint(new p2t::Point(0.7, 0.7), 0);
}

int main() {

    /*
     FIX MEMORY LEAK
     */

    OGRRegisterAll();
    LcpFinder finder{};
    readCostSurface("testpolygon.shp", "targets.shp", "start.shp", &finder);
    //readCostSurfaceDummy("testpolygon.shp", "targets.shp", "start.shp", &finder);
    std::deque<const Coords*> results = finder.leastCostPath();

    std::cout << "LCP SEARCH DONE\n";

    const char *pszDriverName = "ESRI Shapefile";
    OGRSFDriver *poDriver;

    OGRRegisterAll();
    poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
            pszDriverName);
    if (poDriver == NULL) {
        printf("%s driver not available.\n", pszDriverName);
        exit(1);
    }


    OGRDataSource *poDS;
    std::remove("output.shp");
    std::remove("output.shx");
    std::remove("output.prj");
    std::remove("output.dbf");

    poDS = poDriver->CreateDataSource("output.shp", NULL);

    std::cout << "PODS:" << poDS << std::endl;
    if (poDS == NULL) {
        printf("Creation of output file failed.\n");
        exit(1);
    }
    OGRSpatialReference sr;
    sr.importFromEPSG(3047);
    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer("output", &sr, wkbLineString, NULL);
    if (poLayer == NULL) {
        printf("Layer creation failed.\n");
        exit(1);
    }

    OGRFieldDefn oField("cost_end", OFTReal);
    oField.SetPrecision(4);


    if (poLayer->CreateField(&oField) != OGRERR_NONE) {
        printf("Creating Cost field failed.\n");
        exit(1);
    }

    std::map<const Coords*, const Coords*> ancestors;
    std::map<const Coords*, OGRLineString*> geom;


    std::cout << "Initial insert:\n";
    for (const Coords* point : results) {
        ancestors[point] = point;
        geom[point] = new OGRLineString{};
        geom[point]->addPoint(point->getX(), point->getY());
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
            geom[pred] = new OGRLineString{};
            geom[pred]->addPoint(pred->getX(), pred->getY());
            std::cout << "Adding point to pred: " << pred->toString() << std::endl;
        } else {
            ancestors[pred] = ancestor;
            results.push_back(pred);
        }
        geom[ancestor]->addPoint(pred->getX(), pred->getY());
        std::cout << "Adding point to ancestor: " << ancestor->toString() << std::endl;

    }
    std::cout << "Paths mapped\n";



    for (std::pair<const Coords*, OGRLineString*> ls : geom) {
        OGRFeature *poFeature;

        poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
        poFeature->SetField("cost_end", ls.first->getToStart());

        poFeature->SetGeometry(ls.second);

        if (poLayer->CreateFeature(poFeature) != OGRERR_NONE) {
            printf("Failed to create feature in shapefile.\n");
            exit(1);
        }

        OGRFeature::DestroyFeature(poFeature);
        std::cout << "deleting";
        delete ls.second;
    }

    OGRDataSource::DestroyDataSource(poDS);
    return 0;
}
