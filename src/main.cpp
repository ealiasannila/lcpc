//============================================================================
// Name        : lcpc.cpp
// Author      : Elias Annila
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "coords.h"
#include "funnel.h"
#include "min_heap.h"
#include "lcpfinder.h"
#include "geomfunc.h"
#include "../lib/poly2tri.h"
#include <vector>
#include <gdal/gdal.h>
#include <gdal/ogrsf_frmts.h>
#include <ctime>
#include <string.h>
#include<iostream>

bool fileExists(const std::string& name) {
    if (FILE * file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {

        return false;
    }
}

void saveNeighbours(LcpFinder* finder, std::string outfile, Coords f, bool useClosest) {
    const Coords* closest = 0;
    const Coords* c;
    for (auto it = finder->getCoordmap()->begin(); it != finder->getCoordmap()->end(); it++) {
        c = &*it;
        if (closest == 0 or eucDistance(c, &f) < eucDistance(closest, &f)) {
            closest = c;
        }
    }
    nSet n{};
    if (useClosest) {
        for (int pol : closest->belongsToPolygons()) {
            finder->triangulate(pol);
            for (auto it = closest->getNeighbours(pol)->begin(); it != closest->getNeighbours(pol)->end(); it++) {
                std::pair<const Coords*, double> p = *it;
                n.insert(p);
            }
        }
    } else {
        n = finder->findNeighbours(closest);

    }
    OGRSFDriver *driver;
    OGRRegisterAll();
    driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
    OGRDataSource *pointDS;
    if (fileExists(outfile)) {
        driver->DeleteDataSource((outfile).c_str());
    }
    pointDS = driver->CreateDataSource(outfile.c_str(), NULL);
    if (pointDS == NULL) {
        printf("Creation of output file failed.\n");
        exit(1);
    }
    OGRLayer *layer;
    OGRSpatialReference sr;
    sr.importFromEPSG(3047);
    layer = pointDS->CreateLayer("polygon", &sr, wkbLineString, NULL);
    if (layer == NULL) {
        printf("Point layer creation failed.\n");
        exit(1);
    }

    OGRFieldDefn pField("polygon", OFTInteger);


    if (layer->CreateField(&pField) != OGRERR_NONE) {
        printf("Creating polygon field failed.\n");
        exit(1);
    }
    for (std::pair<const Coords*, int> p : n) {
        OGRFeature *feature;
        feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
        OGRLineString line{};
        line.addPoint(closest->getX(), closest->getY());
        line.addPoint(p.first->getX(), p.first->getY());
        feature->SetField("polygon", p.second);
        feature->SetGeometry(&line);
        if (layer->CreateFeature(feature) != OGRERR_NONE) {
            printf("Failed to create feature in shapefile.\n");
            exit(1);
        }
        OGRFeature::DestroyFeature(feature);
    }
    OGRDataSource::DestroyDataSource(pointDS);
    std::cout << "done saving\n";
}

void savePolygon(std::vector<std::vector<std::vector<p2t::Point*>>> polygons, std::string file) {
    if (polygons.empty()) {
        std::cout << "NO POLYGONS";
        exit(1);
    }
    OGRSFDriver *driver;
    OGRRegisterAll();
    driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
    OGRDataSource *pointDS;

    if (fileExists(file)) {
        driver->DeleteDataSource((file).c_str());
    }
    pointDS = driver->CreateDataSource(file.c_str(), NULL);
    if (pointDS == NULL) {
        printf("Creation of output file failed.\n");
        exit(1);
    }
    OGRLayer *pointLayer;
    OGRSpatialReference sr;
    sr.importFromEPSG(3047);
    pointLayer = pointDS->CreateLayer("polygon", &sr, wkbPolygon, NULL);

    if (pointLayer == NULL) {
        printf("Point layer creation failed.\n");
        exit(1);
    }


    for (std::vector<std::vector < p2t::Point*>> polygon : polygons) {
        if (polygon.empty()) {
            continue;
        }
        OGRFeature *poFeature;
        poFeature = OGRFeature::CreateFeature(pointLayer->GetLayerDefn());

        OGRPolygon poly;

        bool ext = true;
        for (std::vector<p2t::Point*> ring : polygon) {
            if (ring.empty()) {
                break;
            }
            OGRLinearRing lr;

            for (p2t::Point* pt : ring) {
                if (pt == 0) {
                    std::cout << "pt ==0\n";
                    exit(1);
                }

                lr.addPoint(pt->x, pt->y);
            }
            lr.addPoint(ring[0]->x, ring[0]->y);
            if (!lr.isClockwise() and ext) {
                std::cout << "reversing\n";
                lr.reverseWindingOrder();

            }
            ext = false;

            poly.addRing(&lr);
        }


        poFeature->SetGeometry(&poly);
        if (pointLayer->CreateFeature(poFeature) != OGRERR_NONE) {
            printf("Failed to create feature in shapefile.\n");
            exit(1);
        }
        OGRFeature::DestroyFeature(poFeature);
    }
    OGRDataSource::DestroyDataSource(pointDS);
    std::cout << "done saving\n";

}

void savePolygons(LcpFinder* finder, std::string outputpolygon) {
    OGRSFDriver *driver;
    OGRRegisterAll();
    driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
    OGRDataSource *pointDS;


    std::cout << outputpolygon << std::endl;
    pointDS = driver->CreateDataSource(outputpolygon.c_str(), NULL);
    if (pointDS == NULL) {
        printf("Creation of output file failed.\n");
        exit(1);
    }
    OGRLayer *pointLayer;
    OGRSpatialReference sr;
    sr.importFromEPSG(3047);
    pointLayer = pointDS->CreateLayer("polygon", &sr, wkbPolygon, NULL);

    if (pointLayer == NULL) {
        printf("Point layer creation failed.\n");
        exit(1);
    }


    for (int i = 0; i < finder->getPolygonCount(); i++) {
        std::vector<std::vector < p2t::Point*>> p2tp = finder->getPolygon(i);

        OGRFeature *poFeature;
        poFeature = OGRFeature::CreateFeature(pointLayer->GetLayerDefn());

        OGRPolygon poly;

        bool ext = true;
        for (std::vector<p2t::Point*> ring : p2tp) {
            OGRLinearRing lr;

            for (p2t::Point* pt : ring) {
                if (pt == 0) {
                    std::cout << "pt ==0\n";
                    exit(1);
                }

                lr.addPoint(pt->x, pt->y);
            }
            lr.addPoint(ring[0]->x, ring[0]->y);
            if (!lr.isClockwise() and ext) {
                std::cout << "reversing\n";
                lr.reverseWindingOrder();

            }
            ext = false;

            poly.addRing(&lr);
        }


        poFeature->SetGeometry(&poly);
        if (pointLayer->CreateFeature(poFeature) != OGRERR_NONE) {
            printf("Failed to create feature in shapefile.\n");
            exit(1);
        }
        OGRFeature::DestroyFeature(poFeature);
    }
    OGRDataSource::DestroyDataSource(pointDS);
    std::cout << "done saving\n";
}

bool comparePolys(OGRPolygon* ogr, std::vector<std::vector<p2t::Point*>> p2tp) {
    if (ogr->getExteriorRing()->getNumPoints() != p2tp[0].size()) {
        std::cout << "ext ring size: " << ogr->getExteriorRing()->getNumPoints() << " vs. " << p2tp[0].size() << std::endl;
        return false;
    }

    if (ogr->getNumInteriorRings() != p2tp.size() - 1) {
        std::cout << "num of int rings: " << ogr->getNumInteriorRings() << " vs. " << p2tp.size() - 1 << std::endl;
        return false;
    }
    for (int i = 0; i < p2tp.size(); i++) {
        if (ogr->getInteriorRing(i)->getNumPoints() != p2tp[i + 1].size()) {
            std::cout << "int ring #" << i + 1 << " ring size: " << ogr->getInteriorRing(i)->getNumPoints() << " vs. " << p2tp[i + 1].size() << std::endl;
            return false;
        }
    }
    return true;

}

void checkCRS(OGRLayer* csLr, OGRLayer* targetLr, OGRLayer* startLr) {
    if (!csLr->GetSpatialRef()->IsProjected()) {
        std::cout << "All source files must be in same projected coordinate system. Geocentric coordinate systems are not supported.\n";
        exit(1);
    }

    int csEpsg = csLr->GetSpatialRef()->GetEPSGGeogCS();
    int targetEpsg = csLr->GetSpatialRef()->GetEPSGGeogCS();
    int startEpsg = csLr->GetSpatialRef()->GetEPSGGeogCS();

    if (csEpsg == -1) {
        std::cout << "WARNING: Unknown cost surface spatial reference EPSG.\n";
    }
    if (targetEpsg == -1) {
        std::cout << "WARNING: Unknown target layer spatial reference EPSG.\n";
    }
    if (startEpsg == -1) {
        std::cout << "WARNING: Unknown start layer spatial reference EPSG.\n";
    }
    if (csEpsg != targetEpsg or csEpsg != startEpsg) {

        std::cout << "ERROR: Cost surface CRS EPSG number doesn't match that of start or target layers. All source files must be in same EPSG projection\n";
        exit(1);
    }


}

void readCostSurface(const char* costSurface, const char* targets, const char* start, LcpFinder* finder, const char* frictionField, std::string maxDist, OGRSpatialReference* sr) {
    double maxd = atof(maxDist.c_str());
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


    checkCRS(csLr, targetLr, startLr);
    *sr = *(csLr->GetSpatialRef());

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
            for (unsigned int ri = 0; ri < csPolygon->getNumInteriorRings(); ri++) {
                OGRLinearRing* intRing = csPolygon->getInteriorRing(ri);

                if (intRing->isClockwise()) {
                    intRing->reverseWindingOrder();
                }
            }
            std::vector<std::vector<std::vector<p2t::Point*> > > sPolygons = simplify(csPolygon);
            for (std::vector<std::vector<p2t::Point*> > polygon : sPolygons) {
                //std::cout<<"POLYGON: "<<pIdx<<"size: "<<polygon[0].size()<<std::endl;
                if (maxd > 0) {
                    intermidiatePoints(&polygon, maxd);
                }
                finder->addPolygon(polygon, csFtre->GetFieldAsDouble(frictionField));
                if (inside(polygon, startp2t)) {
                    finder->addStartPoint(startp2t, pIdx);
                }
                for (int i = 0; i < targetOGRPoints.size(); i++) {
                    p2t::Point* targetp2t = new p2t::Point(targetOGRPoints[i]->getX(), targetOGRPoints[i]->getY());
                    if (inside(polygon, targetp2t)) {
                        finder->addSteinerPoint(targetp2t, pIdx);
                    } else {
                        delete targetp2t;
                    }
                }
                pIdx++;
            }
        } else {
            std::cout << "no polygon geometry " << pIdx << "\n";
            std::cout << csGeometry << std::endl;
            if (csGeometry != 0) {
                std::cout << csGeometry->getGeometryType() << std::endl;
            }
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

void readLinear(const char* linear, LcpFinder* finder, const char* FFFW, const char* FFBW, std::string maxDist) {

    double maxd = atof(maxDist.c_str());
    OGRDataSource *linearDS;
    linearDS = OGRSFDriverRegistrar::Open(linear);
    if (linearDS == NULL) {
        std::cout << linear;
        std::cout << "LINEAR FEATURES NOT FOUND!\n";
        return;
    }
    OGRLayer* linearLr = linearDS->GetLayer(0);
    std::cout << linearLr->GetFeatureCount() << "linear features found" << std::endl;

    OGRFeature * linearFtre;
    while ((linearFtre = linearLr->GetNextFeature()) != NULL) {
        OGRGeometry* linearGeom = linearFtre->GetGeometryRef();
        if (linearGeom != NULL
                && wkbFlatten(linearGeom->getGeometryType()) == wkbLineString) {
            OGRLineString *ls = (OGRLineString *) linearGeom;

            //a point between two first points of linestring. Used to test in which polygon(s) the line string is.
            p2t::Point testPoint{ls->getX(0)+(ls->getX(0) - ls->getX(1)) / 2, ls->getY(0)+(ls->getY(0) - ls->getY(1)) / 2};
            std::array<int, 2> containingPolygons = finder->containingPolygon(&testPoint);
            std::vector<p2t::Point*> linep2t{};
            for (int i = 0; i < ls->getNumPoints(); i++) {
                linep2t.push_back(new p2t::Point{ls->getX(i), ls->getY(i)});
            }
            if (maxd > 0) {
                intermidiatePointsForRing(&linep2t, maxd, false);
            }
            finder->addLine(&linep2t, linearFtre->GetFieldAsDouble(FFFW), linearFtre->GetFieldAsDouble(FFBW), containingPolygons);

            for (p2t::Point* p : linep2t) {
                delete p;
            }

        } else {
            std::cout << "no line geometry " << "\n";
            std::cout << linearGeom << std::endl;
            if (linearGeom != 0) {
                std::cout << linearGeom->getGeometryType() << std::endl;
            }
        }
        OGRFeature::DestroyFeature(linearFtre);
    }
    OGRDataSource::DestroyDataSource(linearDS);
}

bool testDriver(OGRSFDriver* driver) {
    if (driver == NULL) {
        printf("Driver not available Do you want to use default driver ESRI Shapefile?.\n");
        std::string useShp;
        std::cin >> useShp;
        bool isYes = false;
        std::vector<std::string> yesOptions = {"yes", "Yes", "YES", "y", "Y"};
        for (std::string yes : yesOptions) {
            if (useShp == yes) {

                return false;
                break;
            }
        }
        exit(1);
    }
    return true;
}

std::string validateOutfile(OGRSFDriver* driver, std::string outputfile) {

    if (fileExists(outputfile)) {
        std::cout << "Filename already in use, do you want to overwrite? (yes/Yes/YES/Y/y for yes, anything else for no)\n";
        std::string overwrite;
        std::cin >> overwrite;
        bool ow = false;
        std::vector<std::string> yesOptions = {"yes", "Yes", "YES", "y", "Y"};
        for (std::string yes : yesOptions) {
            if (overwrite == yes) {
                driver->DeleteDataSource((outputfile).c_str());
                ow = true;
                break;
            }
        }
        if (!ow) {

            std::cout << "Give new outputfilename that is not in use: \n";
            std::cin >> outputfile;
        }
    }
    return outputfile;
}

void writePoints(std::deque<const Coords*> results, std::string outputfile, std::string pszDriverName, OGRSpatialReference sr, bool overwrite) {
    OGRSFDriver *driver;
    OGRRegisterAll();
    driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName.c_str());

    if (!overwrite) {
        if (!testDriver(driver)) {
            driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
        }
    }


    OGRDataSource *pointDS;

    if (!overwrite) {
        outputfile = validateOutfile(driver, outputfile);
    } else if (fileExists(outputfile)) {
        driver->DeleteDataSource((outputfile).c_str());
    }
    pointDS = driver->CreateDataSource((outputfile).c_str(), NULL);
    if (pointDS == NULL) {
        printf("Creation of output file failed.\n");
        exit(1);
    }


    OGRLayer *pointLayer;
    pointLayer = pointDS->CreateLayer("point_results", &sr, wkbPoint, NULL);
    if (pointLayer == NULL) {
        printf("Point layer creation failed.\n");
        exit(1);
    }

    OGRFieldDefn cField("cost", OFTReal);
    cField.SetPrecision(2);


    if (pointLayer->CreateField(&cField) != OGRERR_NONE) {
        printf("Creating Cost field failed.\n");
        exit(1);
    }

    for (const Coords* point : results) {
        OGRFeature *poFeature;

        poFeature = OGRFeature::CreateFeature(pointLayer->GetLayerDefn());
        poFeature->SetField("cost", point->getToStart());
        OGRPoint ogrpt(point->getX(), point->getY());
        poFeature->SetGeometry(&ogrpt);

        if (pointLayer->CreateFeature(poFeature) != OGRERR_NONE) {

            printf("Failed to create feature in shapefile.\n");
            exit(1);
        }

        OGRFeature::DestroyFeature(poFeature);
    }

    OGRDataSource::DestroyDataSource(pointDS);


}

void writeOutput(std::deque<const Coords*> results, std::string outputfile, std::string pszDriverName, OGRSpatialReference sr, bool overwrite) {

    OGRSFDriver *driver;
    OGRRegisterAll();
    driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName.c_str());
    if (!overwrite) {
        if (!testDriver(driver)) {
            driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
        }
    }
    OGRDataSource *pathDS;
    if (!overwrite) {
        outputfile = validateOutfile(driver, outputfile);

    } else if (fileExists(outputfile)) {
        driver->DeleteDataSource((outputfile).c_str());
    }
    pathDS = driver->CreateDataSource((outputfile).c_str(), NULL);
    if (pathDS == NULL) {
        printf("Creation of output file failed.\n");
        exit(1);
    }


    OGRLayer *pathLayer;
    pathLayer = pathDS->CreateLayer("least_cost_path", &sr, wkbLineString, NULL);
    if (pathLayer == NULL) {
        printf("Path layer creation failed.\n");
        exit(1);
    }

    OGRFieldDefn cField("cost_end", OFTReal);
    cField.SetPrecision(2);


    OGRFieldDefn xField("x", OFTReal);
    xField.SetPrecision(2);

    OGRFieldDefn yField("y", OFTReal);
    yField.SetPrecision(2);


    if (pathLayer->CreateField(&cField) != OGRERR_NONE) {
        printf("Creating Cost field failed.\n");
        exit(1);
    }
    if (pathLayer->CreateField(&xField) != OGRERR_NONE) {
        printf("Creating Cost field failed.\n");
        exit(1);
    }
    if (pathLayer->CreateField(&yField) != OGRERR_NONE) {

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



    for (const Coords* point : results) {



        minheap.push(point);
        ancestors[point] = point;
        geom[point] = new OGRLineString{};
        geom[point]->addPoint(point->getX(), point->getY());
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

    for (std::pair<const Coords*, OGRLineString*> ls : geom) {
        OGRFeature *poFeature;

        poFeature = OGRFeature::CreateFeature(pathLayer->GetLayerDefn());
        poFeature->SetField("cost_end", ls.first->getToStart());
        poFeature->SetField("x", ls.first->getX());
        poFeature->SetField("y", ls.first->getY());


        poFeature->SetGeometry(ls.second);

        if (pathLayer->CreateFeature(poFeature) != OGRERR_NONE) {

            printf("Failed to create feature in shapefile.\n");
            exit(1);
        }

        OGRFeature::DestroyFeature(poFeature);
        delete ls.second;
    }

    OGRDataSource::DestroyDataSource(pathDS);

}

void writeOutputInvidual(std::deque<const Coords*> results, std::string outputfile, std::string pszDriverName, OGRSpatialReference sr, bool overwrite) {

    OGRSFDriver *driver;
    OGRRegisterAll();
    driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName.c_str());
    if (!overwrite) {
        if (!testDriver(driver)) {
            driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
        }
    }
    OGRDataSource *pathDS;
    if (!overwrite) {
        outputfile = validateOutfile(driver, outputfile);

    } else if (fileExists(outputfile)) {
        driver->DeleteDataSource((outputfile).c_str());
    }
    pathDS = driver->CreateDataSource((outputfile).c_str(), NULL);
    if (pathDS == NULL) {
        printf("Creation of output file failed.\n");
        exit(1);
    }


    OGRLayer *pathLayer;
    pathLayer = pathDS->CreateLayer("least_cost_path", &sr, wkbLineString, NULL);
    if (pathLayer == NULL) {
        printf("Path layer creation failed.\n");
        exit(1);
    }

    OGRFieldDefn cField("cost", OFTReal);
    cField.SetPrecision(2);

    OGRFieldDefn idField("path_id", OFTInteger);


    if (pathLayer->CreateField(&cField) != OGRERR_NONE) {
        printf("Creating Cost field failed.\n");
        exit(1);
    }
    if (pathLayer->CreateField(&idField) != OGRERR_NONE) {
        printf("Creating Cost field failed.\n");
        exit(1);
    }

    int i = 1;
    for (const Coords* point : results) {
        OGRFeature *poFeature;

        poFeature = OGRFeature::CreateFeature(pathLayer->GetLayerDefn());
        poFeature->SetField("cost", point->getToStart());
        poFeature->SetField("path_id", i);
        i++;

        OGRLineString ogrls{};
        const Coords* next = point;
        while (next != 0) {
            ogrls.addPoint(next->getX(), next->getY());
            next = next->getPred();
        }
        poFeature->SetGeometry(&ogrls);

        if (pathLayer->CreateFeature(poFeature) != OGRERR_NONE) {
            printf("Failed to create feature in shapefile.\n");
            exit(1);
        }

        OGRFeature::DestroyFeature(poFeature);

    }
    OGRDataSource::DestroyDataSource(pathDS);
}

bool argExists(std::string argSwitch, char* argv[], int argc) {
    for (int i = 0; i < argc; i++) {
        if (argSwitch.compare(argv[i]) == 0) {
            return true;
        }
    }
    return false;
}

std::string getArgVal(std::string argSwitch, char* argv[], int argc) {
    for (int i = 0; i < argc; i++) {
        if (argSwitch.compare(argv[i]) == 0) {
            return argv[i + 1];
        }
    }
}

void printFinder(LcpFinder* finder) {
    std::cout << finder->getCoordmap()->size() << std::endl;
}

struct PointDistComparator {
    p2t::Point* p;

    PointDistComparator(p2t::Point* p) {
        this->p = p;
    }

    inline bool operator()(p2t::Point* pt1, p2t::Point* pt2) {
        return (eucDistance(p, pt1) < eucDistance(p, pt2));
    }
};

int checkIntersections(std::vector<OGRLineString*> lines, p2t::Point* p, int ringI, OGRLinearRing* ring, int ringSize, int nextI, std::vector<p2t::Point*>* extRing, std::deque < std::array<int, 2 >> *vertexQue, std::vector<std::vector < p2t::Point*>>*extras) {

    OGRLineString* closestLS = 0;
    double closestD = std::numeric_limits<double>::max();
    bool fromEntry = true;

    p2t::Point np{ring->getX(nextI), ring->getY(nextI)};
    for (OGRLineString* ls : lines) {
        int lineSize = ls->getNumPoints();
        p2t::Point entry{ls->getX(0), ls->getY(0)};
        p2t::Point exit{ls->getX(lineSize - 1), ls->getY(lineSize - 1)};
        double entryd = eucDistance(p, &entry);
        double exitd = eucDistance(p, &exit);

        if (pointOnSegment(p, &np, &entry)and entryd < closestD and entryd > 0) {
            closestD = entryd;
            fromEntry = true;
            closestLS = ls;
        }
        if (pointOnSegment(p, &np, &exit) and exitd < closestD and exitd > 0) {
            closestD = exitd;
            fromEntry = false;
            closestLS = ls;
        }
    }
    if (closestLS != 0) {
        vertexQue->push_back(std::array<int, 2>{ringI, nextI});
        int start;
        int end;
        int step;
        if (fromEntry) {
            start = 0;
            end = closestLS->getNumPoints();
            step = 1;
        } else {
            start = closestLS->getNumPoints() - 1;
            step = -1;
            end = -1;
        }
        for (int li = start; li != end; li += step) {
            extRing->push_back(new p2t::Point{closestLS->getX(li), closestLS->getY(li)});
        }

        int newNext = (previousPoint(ring, extRing->back()) + 1) % ringSize;
        
        //intersecting line doubles back to same edge.
        if (newNext == nextI) {
            std::vector<p2t::Point*> newExterior{};
            for (int li = end - step; li != start - step; li -= step) {
                newExterior.push_back(new p2t::Point{closestLS->getX(li), closestLS->getY(li)});
            }
            extras->push_back(newExterior);
            std::cout<<"same edge\n";
        }
        //intersecting lines cross two same edges
        
        return newNext;
    }
    return -1;
}

void preprocess(const char* costSurface, const char* linear) {
    OGRDataSource *linearDS;
    linearDS = OGRSFDriverRegistrar::Open(linear);
    if (linearDS == NULL) {
        std::cout << linear;
        std::cout << " LINEAR FEATURES NOT FOUND!\n";
        exit(1);
    }
    OGRLayer* linearLr = linearDS->GetLayer(0);
    std::cout << linearLr->GetFeatureCount() << " linear features found" << std::endl;

    OGRDataSource *csDS;
    csDS = OGRSFDriverRegistrar::Open(costSurface);
    if (csDS == NULL) {
        std::cout << costSurface;
        std::cout << " COST SURFACENOT FOUND!\n";
        exit(1);
    }
    OGRLayer* csLr = csDS->GetLayer(0);
    std::cout << csLr->GetFeatureCount() << " cost areas found" << std::endl;

    OGRFeature * csFtre;
    std::vector<std::vector<std::vector < p2t::Point*>>> outPolygons;
    while ((csFtre = csLr->GetNextFeature()) != NULL) {
        OGRGeometry* csGeom = csFtre->GetGeometryRef();
        if (csGeom != NULL && wkbFlatten(csGeom->getGeometryType()) == wkbPolygon) {
            OGRPolygon *polygon = (OGRPolygon *) csGeom;


            //ADD CASE IF CONTAINS
            linearLr->ResetReading();
            OGRFeature * linearFtre;
            std::vector<OGRLineString*> lines{};
            while ((linearFtre = linearLr->GetNextFeature()) != NULL) {
                OGRGeometry* linearGeom = linearFtre->GetGeometryRef();
                if (linearGeom != NULL && wkbFlatten(linearGeom->getGeometryType()) == wkbLineString) {
                    OGRLineString *ls = (OGRLineString *) linearGeom;
                    if (polygon->Intersects(ls)) {
                        OGRGeometry* intersection = polygon->Intersection(ls);
                        if (wkbFlatten(intersection->getGeometryType()) == wkbLineString) {
                            lines.push_back((OGRLineString*) intersection);
                        } else if (wkbFlatten(intersection->getGeometryType()) == wkbMultiLineString) {
                            OGRMultiLineString* multiline = (OGRMultiLineString*) intersection;
                            for (int i = 0; i < multiline->getNumGeometries(); i++) {
                                lines.push_back((OGRLineString*) multiline->getGeometryRef(i));
                            }
                        }
                    }

                }
            }
            std::cout << std::fixed;
            //line can intersect at most start and endpoint. Must ińtersect at either
            if (!lines.empty()) {
                std::deque < std::array<int, 2 >> vertexQue
                {
                };
                vertexQue.push_back(std::array<int, 2>{0, 0});

                std::set<std::array<int, 2 >> used;
                while (!vertexQue.empty()) {
                    int i = vertexQue.front()[1];
                    int ringI = vertexQue.front()[0];
                    vertexQue.pop_front();
                    std::vector<p2t::Point*> extRing{};
                    int limit = 0;
                    while (limit < 100 and used.find(std::array<int, 2>{ringI, i}) == used.end()) {
                        used.insert(std::array<int, 2>{ringI, i});

                        limit++;
                        OGRLinearRing* ring;

                        if (ringI == 0) {
                            ring = polygon->getExteriorRing();
                        } else {
                            ring = polygon->getInteriorRing(ringI - 1);
                        }
                        int ringSize = ring->getNumPoints();
                        if (abs(ring->getX(0) - ring->getX(ringSize - 1)) < 0.00001 and abs(ring->getY(0) - ring->getY(ringSize - 1) < 0.00001)) {
                            ringSize--;
                        }
                        int nextI = (i + 1) % ringSize;


                        p2t::Point* p = new p2t::Point{ring->getX(i), ring->getY(i)};
                        extRing.push_back(p);
                        std::vector<std::vector < p2t::Point*>> extras;
                        int ci = -1;
                        do {
                            ci = checkIntersections(lines, extRing.back(), ringI, ring, ringSize, nextI, &extRing, &vertexQue, &extras);
                            if (ci != -1) {
                                nextI = ci;
                            }
                        } while (ci != -1);
                        i = nextI;
                        std::cout << extras.size();
                        for (std::vector<p2t::Point*> extra : extras) {

                            outPolygons.push_back(std::vector<std::vector < p2t::Point*>>
                            {
                                extra
                            });
                        }

                    }
                    if (!extRing.empty()) {
                        outPolygons.push_back(std::vector<std::vector < p2t::Point*>>
                        {
                            extRing
                        });
                    }
                }
            }
        }
    }
    savePolygon(outPolygons, "testdata/out.shp");

}

int main(int argc, char* argv[]) {

    OGRRegisterAll();
    preprocess("testdata/test_sim.shp", "testdata/linear.shp");
    exit(0);
    if (argc < 2 or strcmp(argv[1], "-h") == 0) {
        std::cout << "This program is used to search for least cost paths in polygonal costsurface for more information see: URL.\n";
        std::cout << "USAGE:\n"
                "The program requires at least 4 parameters to run correctly. Parameters should be in following order:\n"
                "\tname of cost surface file. (Supports formats supported by OGR)\n"
                "\tname of target points file can be any number of points. All points must be inside cost surface polygons. (Supports formats supported by OGR)\n"
                "\tname of start point file, must be single point inside one of the cost surface polygons. (Supports formats supported by OGR)\n"
                "\tname of the friction field in cost surface\n"
                "Addittionally following parameters can be specified with formt <switch> <value> (for example -o out_path.shp):\n"
                "\t -o output_path_file_name, if no extension is given a folder will be assumed (at least with shapefile driver)\n"
                "\t -i output_path_file_name, returns outputpaths as invidual paths with duplicates."
                "\t -p output_points_file_name, if no extension is given a folder will be assumed (at least with shapefile driver)\n"
                "\t -d maximum_distance_between_nodes in cost surface. Default value is 0 (no nodes added). This is used to add temporary additional nodes during LCP calculation if nodes are too far apart.\n"
                "\t -l linear cost features\n"
                "\t --fwc forwards cost field\n"
                "\t --bwc backwards cost field\n"
                "\t --driver output_driver_name. Default value is ESRI Shapefile (Supports drivers supported by OGR)\n"
                "\t -a algorithm to use in lcp search. Astar is default value. Options are dijksra/astar. Using astar is recommended unless number of targetpoints is high compared to number of nodes."
                "\t --overwrite if used will not prompt on anything and just overwrite all files";
        exit(0);
    }

    /*
    std::cout << argc << " provided:";
    for (int i = 1; i < argc; i++) {
        std::cout << (i) << ": " << argv[i] << std::endl;
    }
     */

    if (argc < 5) {
        std::cout << "Invalid arguments provided see \"lcp -h\" for details\n";
        return 0;
    }
    int alg = 1;
    if (argExists("-a", argv, argc)) {
        std::string algArg = getArgVal("-a", argv, argc);
        if (algArg.compare("astar") == 0) {
            alg = 1;
        } else if (algArg.compare("dijkstra") == 0) {
            alg = 0;
        } else {
            std::cout << "unknown search algorithm. using A* (possible choises are \"astar\" or \"dijkstra\"\n";
        }
    }



    bool overwrite = argExists("--overwrite", argv, argc);
    std::string distance = "0";
    if (argExists("-d", argv, argc)) {
        distance = getArgVal("-d", argv, argc);
    }

    LcpFinder finder{};
    std::cout << "Reading cost surface...\n";
    OGRSpatialReference sr;
    std::clock_t begin = std::clock();

    readCostSurface(argv[1], argv[2], argv[3], &finder, argv[4], distance, &sr);
    double secs = double(std::clock() - begin) / CLOCKS_PER_SEC;
    if (argExists("--scs", argv, argc)) {
        std::string outcs = getArgVal("--scs", argv, argc);
        std::cout << "Saving modified costsurface\n";
        savePolygons(&finder, outcs);
    }
    //std::cout << "saving neighburs" << std::endl;
    //saveNeighbours(&finder, "testdata/neighbours.shp", Coords{337018, 6704919}, false);
    //exit(0);
    if (argExists("-l", argv, argc)) {
        readLinear(getArgVal("-l", argv, argc).c_str(), &finder, getArgVal("--fwc", argv, argc).c_str(), getArgVal("--bwc", argv, argc).c_str(), distance);
    }
    std::cout << "Finished reading cost surface (took " << secs << " s). Starting LCP search...\n";
    begin = std::clock();
    std::deque<const Coords*> results = finder.leastCostPath(alg);
    secs = double(std::clock() - begin) / CLOCKS_PER_SEC;
    std::cout << "\nSearch finished (took " << secs << "s). Writing results...\n";

    std::string driver = "ESRI Shapefile";
    if (argExists("--driver", argv, argc)) {
        driver = getArgVal("--driver", argv, argc);
    }

    if (argExists("-o", argv, argc)) {
        writeOutput(results, getArgVal("-o", argv, argc), driver, sr, overwrite);
    }

    if (argExists("-i", argv, argc)) {
        writeOutputInvidual(results, getArgVal("-i", argv, argc), driver, sr, overwrite);
    }
    if (argExists("-p", argv, argc)) {
        writePoints(results, getArgVal("-p", argv, argc), driver, sr, overwrite);
    }
    std::cout << "All done!\n";
    return 0;
}
