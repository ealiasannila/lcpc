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
#include <ctime>
#include <string.h>
#include <stdlib.h>

bool inside(OGRPolygon* polygon, OGRPoint* point) {
    if (polygon->getExteriorRing()->isPointInRing(point)) {
        for (unsigned int ri = 0; ri < polygon->getNumInteriorRings(); ri++) {
            if (polygon->getInteriorRing(ri)->isPointInRing(point)) {
                return false;
            }
        }
        return true;
    }
    return false;
}

void readCostSurface(const char* costSurface, const char* targets, const char* start, LcpFinder* finder, const char* frictionField, const char* maxDist) {
    double maxd = atof(maxDist);
    OGRDataSource *csDS;
    csDS = OGRSFDriverRegistrar::Open(costSurface);
    if (csDS == NULL) {

        std::cout << costSurface;
        std::cout << " COST SURFACE NOT FOUND!\n";
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
                    p2t::Point* point = new p2t::Point(intRing->getX(i), intRing->getY(i));
                    polygon[ri + 1].push_back(point);
                }
                delete polygon[ri + 1].back();
                polygon[ri + 1].pop_back();
            }

            intermidiatePoints(&polygon, maxd);

            finder->addPolygon(polygon, csFtre->GetFieldAsDouble(frictionField));

            if (inside(csPolygon, startPoint)) {
                finder->addStartPoint(startp2t, pIdx);
            }
            for (int i = 0; i < targetOGRPoints.size(); i++) {
                if (inside(csPolygon, targetOGRPoints[i])) {
                    finder->addSteinerPoint(new p2t::Point(targetOGRPoints[i]->getX(), targetOGRPoints[i]->getY()), pIdx);
                }
            }


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

void writeShapeFile(std::deque<const Coords*> results, std::string outputfile) {

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
    std::remove((outputfile + ".shp").c_str());
    std::remove((outputfile + ".shx").c_str());
    std::remove((outputfile + ".prj").c_str());
    std::remove((outputfile + ".dbf").c_str());

    poDS = poDriver->CreateDataSource((outputfile + ".shp").c_str(), NULL);

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

    OGRFieldDefn cField("cost_end", OFTReal);
    cField.SetPrecision(2);


    OGRFieldDefn xField("x", OFTReal);
    xField.SetPrecision(2);

    OGRFieldDefn yField("y", OFTReal);
    yField.SetPrecision(2);


    if (poLayer->CreateField(&cField) != OGRERR_NONE) {
        printf("Creating Cost field failed.\n");
        exit(1);
    }
    if (poLayer->CreateField(&xField) != OGRERR_NONE) {
        printf("Creating Cost field failed.\n");
        exit(1);
    }
    if (poLayer->CreateField(&yField) != OGRERR_NONE) {
        printf("Creating Cost field failed.\n");
        exit(1);
    }

    std::map<const Coords*, const Coords*> ancestors;
    std::map<const Coords*, OGRLineString*> geom;

    struct compToStart {

        bool operator()(const Coords* x, const Coords* y) const {
            return ((x->getToStart()) < (y->getToStart()));
        }
    } compare{};
    MinHeap<const Coords*, compToStart> minheap(compare);



    std::cout << "Initial insert:\n";
    for (const Coords* point : results) {
        minheap.push(point);
        ancestors[point] = point;
        geom[point] = new OGRLineString{};
        geom[point]->addPoint(point->getX(), point->getY());
        std::cout << point->toString() << std::endl;
    }

    while (!minheap.empty()) {
        const Coords* point = minheap.top();
        const Coords* pred = point->getPred();
        minheap.pop();
        if (pred == 0) {
            continue;
        }

        const Coords* ancestor = ancestors[point];

        if (ancestors.find(pred) != ancestors.end() and pred->getPred() != 0) {
            ancestors[pred] = pred;
            geom[pred] = new OGRLineString{};
            geom[pred]->addPoint(pred->getX(), pred->getY());
        } else {
            ancestors[pred] = ancestor;
            minheap.push(pred);
        }
        geom[ancestor]->addPoint(pred->getX(), pred->getY());

    }
    std::cout << "Paths mapped\n";



    for (std::pair<const Coords*, OGRLineString*> ls : geom) {
        OGRFeature *poFeature;

        poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
        poFeature->SetField("cost_end", ls.first->getToStart());
        poFeature->SetField("x", ls.first->getX());
        poFeature->SetField("y", ls.first->getY());


        poFeature->SetGeometry(ls.second);

        if (poLayer->CreateFeature(poFeature) != OGRERR_NONE) {
            printf("Failed to create feature in shapefile.\n");
            exit(1);
        }

        OGRFeature::DestroyFeature(poFeature);
        delete ls.second;
    }

    OGRDataSource::DestroyDataSource(poDS);
}

int main(int argc, char* argv[]) {
    if (argc != 7) {
        std::cout << "Invalid arguments provided\n";
        std::cout << "Correct form is lcpc cost_surface.shp target.shp start.shp name_of_friction_field output_file_name_without_shp_extension max_distance_between_nodes\n";
        return 0;
    }

    OGRRegisterAll();
    LcpFinder finder{};
    readCostSurface(argv[1], argv[2], argv[3], &finder, argv[4], argv[6]);
    //readCostSurfaceDummy("testpolygon.shp", "targets.shp", "start.shp", &finder);
    std::cout << "READ SURFACE DONE\n";

    std::deque<const Coords*> results = finder.leastCostPath();

    std::cout << "TIME SPENT IN FUNNELALGORITHM: " << finder.funnel_secs << " SECONDS\n";
    std::cout << "-funnelqueing: " << finder.fq_secs << " SECONDS\n";
    std::cout << "-getting opposing: " << finder.base_secs << " SECONDS\n";
    std::cout << "-reacting: " << finder.react_secs << " SECONDS\n";

    writeShapeFile(results, argv[5]);
    std::cout << "All done!\n";
    return 0;
}
