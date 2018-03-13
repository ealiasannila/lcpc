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
#include <sstream>
#include <map>

bool fileExists(const std::string& name) {
    if (FILE * file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {

        return false;
    }
}

void saveTriangulation(LcpFinder* finder, std::string outfile, int polygon) {
    std::set<const Triangle*> used;
    finder->triangulate(polygon);

    GDALDriver *driver;
    GDALAllRegister();
    driver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
    GDALDataset *pointDS;
    if (fileExists(outfile)) {
        driver->Delete((outfile).c_str());
    }
    pointDS = driver->Create(outfile.c_str(), NULL, NULL, NULL, GDT_Unknown, NULL);
    if (pointDS == NULL) {
        printf("Creation of output file failed.\n");
        exit(1);
    }
    OGRLayer *layer;
    OGRSpatialReference sr;
    sr.importFromEPSG(3047);
    layer = pointDS->CreateLayer("polygon", &sr, wkbPolygon, NULL);
    if (layer == NULL) {
        printf("Point layer creation failed.\n");
        exit(1);
    }

    OGRFieldDefn pField("polygon", OFTInteger);


    if (layer->CreateField(&pField) != OGRERR_NONE) {
        printf("Creating polygon field failed.\n");
        exit(1);
    }
    int tricount = 0;
    int nodecount = 0;
    std::vector<std::vector<const Coords*> > pol = finder->getPolygon(polygon);
    for (std::vector<const Coords*> ring : pol) {
        for (auto it = ring.begin(); it != ring.end(); it++) {
            const Coords c = **it;
            nodecount++;
            if (c.belongsToPolygon(polygon) && c.flag == 0) {
                for (auto p = c.getTriangles(polygon)->begin(); p != c.getTriangles(polygon)->end(); p++) {
                    OGRFeature *feature;
                    feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
                    OGRPolygon polygon{};
                    OGRLinearRing ring{};
                    if (used.find(*p) == used.end()) {
                        tricount++;
                        for (int i = 0; i < 3; i++) {
                            ring.addPoint((*p)->points[i]->getX(), (*p)->points[i]->getY());
                        }
                        ring.addPoint((*p)->points[0]->getX(), (*p)->points[0]->getY());
                        used.insert(*p);
                    }
                    polygon.addRing(&ring);
                    feature->SetGeometry(&polygon);
                    if (layer->CreateFeature(feature) != OGRERR_NONE) {
                        printf("Failed to create feature in shapefile.\n");
                        exit(1);
                    }
                    OGRFeature::DestroyFeature(feature);
                }
            }
        }
    }
    std::cout << "counts:" << nodecount << "  " << tricount << std::endl;

    GDALClose(pointDS);
    std::cout << "done saving\n";




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

    n = finder->findNeighbours(closest);


}
 

/*
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
 */

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

void readCostSurface(const char* costSurface, const char* targets, const char* start, LcpFinder* finder, const char* frictionField, OGRSpatialReference* sr) {

    GDALDataset *csDS;
    GDALDataset *startDS;
    GDALDataset *targetDS;

    csDS = (GDALDataset*) GDALOpenEx(costSurface, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (csDS == NULL) {
        printf("Cost surface not found.\n");
        return;
    }

    targetDS = (GDALDataset*) GDALOpenEx(targets, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (targetDS == NULL) {
        printf("Cost surface not found.\n");
        return;
    }

    startDS = (GDALDataset*) GDALOpenEx(start, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (targetDS == NULL) {
        printf("Cost surface not found.\n");
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

    OGRFeature *startFtre = startLr->GetFeature(0);
    std::vector<OGRPoint*> targetOGRPoints;
    std::vector<OGRFeature*> targetPointers;
    OGRFeature * targetFtre;
    while ((targetFtre = targetDS->GetNextFeature(&targetLr, NULL, NULL, NULL)) != NULL) {
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
    while ((csFtre = csDS->GetNextFeature(&csLr, NULL, NULL, NULL)) != NULL) {
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
            std::vector<std::vector<std::vector<p2t::Point> > > sPolygons = simplify(csPolygon);
            for (std::vector<std::vector<p2t::Point> > polygon : sPolygons) {
                //std::cout<<"POLYGON: "<<pIdx<<"size: "<<polygon[0].size()<<std::endl;
                finder->addPolygon(polygon, csFtre->GetFieldAsDouble(frictionField));
                if (inside(polygon, *startp2t)) {
                    finder->addStartPoint(startp2t, pIdx);

                }
                for (int i = 0; i < targetOGRPoints.size(); i++) {
                    p2t::Point* targetp2t = new p2t::Point(targetOGRPoints[i]->getX(), targetOGRPoints[i]->getY());
                    if (inside(polygon, *targetp2t)) {
                        finder->addTargetPoint(targetp2t, pIdx);
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

    GDALClose(csDS);
    GDALClose(startDS);
    GDALClose(targetDS);
}

void readLinear(const char* linear, LcpFinder* finder, const char* FFFW, const char* FFBW, const char* FFC, double maxd) {

    GDALDataset *linearDS;
    linearDS = (GDALDataset*) GDALOpenEx(linear, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (linearDS == NULL) {
        std::cout << linear;
        std::cout << "LINEAR FEATURES NOT FOUND!\n";
        return;
    }
    OGRLayer* linearLr = linearDS->GetLayer(0);
    std::cout << linearLr->GetFeatureCount() << "linear features found" << std::endl;

    OGRFeature * linearFtre;
    while ((linearFtre = linearDS->GetNextFeature(&linearLr, NULL, NULL, NULL)) != NULL) {

        OGRGeometry* linearGeom = linearFtre->GetGeometryRef();
        if (linearGeom != NULL
                && wkbFlatten(linearGeom->getGeometryType()) == wkbLineString) {
            OGRLineString *ls = (OGRLineString *) linearGeom;
            std::vector<p2t::Point*> linep2t{};
            for (int i = 0; i < ls->getNumPoints(); i++) {
                linep2t.push_back(new p2t::Point{ls->getX(i), ls->getY(i)});
            }

            if (maxd > 0) {
                intermidiatePointsForRing(&linep2t, maxd, false);
            }

            finder->addLine(&linep2t, linearFtre->GetFieldAsDouble(FFFW), linearFtre->GetFieldAsDouble(FFC));



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

bool testDriver(GDALDriver* driver) {
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

std::string validateOutfile(GDALDriver* driver, std::string outputfile) {

    if (fileExists(outputfile)) {
        std::cout << "Filename already in use, do you want to overwrite? (yes/Yes/YES/Y/y for yes, anything else for no)\n";
        std::string overwrite;
        std::cin >> overwrite;
        bool ow = false;
        std::vector<std::string> yesOptions = {"yes", "Yes", "YES", "y", "Y"};
        for (std::string yes : yesOptions) {
            if (overwrite == yes) {
                driver->pfnDeleteDataSource(driver, (outputfile).c_str());
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

/*

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
 */
void writeOutput(std::deque<const Coords*> results, std::string outputfile, std::string pszDriverName, OGRSpatialReference sr, bool overwrite) {
    GDALDriver *driver;
    GDALAllRegister();
    driver = GetGDALDriverManager()->GetDriverByName(pszDriverName.c_str());
    if (!overwrite) {
        if (!testDriver(driver)) {
            driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
        }
    }
    GDALDataset *pathDS;
    if (!overwrite) {
        outputfile = validateOutfile(driver, outputfile);

    } else if (fileExists(outputfile)) {
        driver->Delete((outputfile).c_str());
    }
    pathDS = driver->Create((outputfile).c_str(), NULL, NULL, NULL, GDT_Unknown, NULL);
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

    GDALClose(pathDS);
    std::cout<<"wrote output to "<<outputfile<<std::endl;
}

/*
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
 * 
 * */

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

void printTriangle(p2t::Triangle* t) {
    std::cout << "Triangle at " << t << "\n";
    if (t != 0) {
        for (int i = 0; i < 3; i++) {
            std::cout << t->GetPoint(i)->x << "," << t->GetPoint(i)->y << std::endl;
        }
    } else {
        std::cout << "NULL\n";
    }
}

void printTriangle(Triangle* t) {
    std::cout << "Triangle : " << t << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "p: " << t->points[i]->getX() << "," << t->points[i]->getY() << std::endl;
        std::cout << "n: " << t->neighbours[i] << std::endl;
    }
}

int main(int argc, char* argv[]) {

    /*
    std::string ext = "773.6664582066915 -1097.654455639297, 1438.37406931964 -1164.700449293967, 1980.179460847239 -673.1656611039796, 554.494146924882 -772.2186036729893, 1074.571231374806 -428.9340870775081, 2264.825880576197 -257.0309796962655, 3220.187676508343 -349.1995507060333, 3074.961489088575 -729.02188703466, 2728.652888318357 -656.4087933247756, 2533.156097560975 -1287.584146341464, 3879.291142490372 -2533.177984595636, 3300.489063671967 -2982.169162365442, 2421.443645699615 -2784.531001283698, 2695.139152759949 -1941.101989730424, 2582.766867233537 -1802.211427409336, 1438.37406931964 -1667.40648267009, 1202.319938098457 -1964.937687559853, 680.8580456852442 -1885.975691857102, 382.6913992297809 -1326.68350449294";
    std::string inter = "1661.910108269481 -1163.02773944267, 1017.638571286707 -1516.94296409984, 889.7133537179041 -1516.94296409984, 889.7133537179041 -1625.883521049343, 1783.871210543036 -1431.967605994615, 2018.411791838334 -1535.165461764547, 1852.669781056323 -1372.55065873314, 1877.687443061155 -1277.170822339719, 1994.957733708804 -1341.2785812271, 2317.060132021013 -1549.237896642265, 2248.261561507727 -1259.971179711397, 2130.991270860077 -1133.319265811936, 2367.095456030677 -1061.393487548045, 2113.791628231755 -915.9783271449601, 1993.394129833502 -953.5048201522078, 1946.486013574443 -1163.02773944267";

    std::istringstream isext(ext);
    std::istringstream isinter(inter);

    std::vector<std::vector <const Coords*>> polygon(1);
    while (isext.good()) {
        std::string substr;
        std::string::size_type sz; // alias of size_t

        std::getline(isext, substr, ',');

        int i = substr.find(" ");
        double x = std::stod(substr, &sz);
        double y = std::stod(substr.substr(sz));

        polygon[0].push_back(new Coords(x, y, 1));
    }


    while (isinter.good()) {
        std::string substr;
        std::getline(isinter, substr, ',');
        std::string::size_type sz; // alias of size_t

        int i = substr.find(" ");
        double x = std::stod(substr, &sz);
        double y = std::stod(substr.substr(sz));

    //    polygon[1].push_back(new Coords(x, y, 1));
    }
    //std::reverse(polygon[1].begin(), polygon[1].end());
    splitToMonotone(polygon);

    return 0;
     * */
    GDALAllRegister();

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
    double distanceVal = atof(distance.c_str());
    distanceVal = 25;
    finder.setMaxD(distanceVal);
    std::cout << "Reading cost surface...\n";
    OGRSpatialReference sr;
    std::clock_t begin = std::clock();

    readCostSurface(argv[1], argv[2], argv[3], &finder, argv[4], &sr);
    double secs = double(std::clock() - begin) / CLOCKS_PER_SEC;
    if (argExists("--scs", argv, argc)) {
        std::string outcs = getArgVal("--scs", argv, argc);
        std::cout << "Saving modified costsurface\n";
        //      savePolygons(&finder, outcs);
    }
    //std::cout << "saving neighburs" << std::endl;

    if (argExists("-l", argv, argc)) {
        std::cout << "Reading linear\n";
        std::cout << getArgVal("--fwc", argv, argc).c_str()<< std::endl;
        std::cout << getArgVal("--bwc", argv, argc).c_str()<< std::endl;
        std::cout << getArgVal("--cc", argv, argc).c_str()<< std::endl;

        readLinear(getArgVal("-l", argv, argc).c_str(), &finder, getArgVal("--fwc", argv, argc).c_str(), getArgVal("--bwc", argv, argc).c_str(), getArgVal("--cc", argv, argc).c_str(), distanceVal);
        //finder.addBuffers(1);
    }
    std::cout << "Finished reading cost surface (took " << secs << " s). Starting LCP search...\n";
    //saveNeighbours(&finder, "testdata/neighbours_2.shp", Coords{305865.9,6719385.9}, false);
    //exit(0);
    /*
    for (int i = 0; i < finder.getPolygonCount(); i++) {
        std::cout << "polygon: " << i << "/" << finder.getPolygonCount() << std::endl;
        finder.triangulate(i);
    }*/
    //saveTriangulation(&finder, "testdata/triangulation.shp", 12);
    //finder.triangulate(251);
    //exit(0);
    begin = std::clock();
    std::deque<const Coords*> results = finder.leastCostPath(alg);
    //std::cout<<"found: "<<results.size()<<std::endl;
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
        //writeOutputInvidual(results, getArgVal("-i", argv, argc), driver, sr, overwrite);
    }

    if (argExists("-p", argv, argc)) {
        //writePoints(results, getArgVal("-p", argv, argc), driver, sr, overwrite);
    }
    std::cout << "All done!\n";

    return 0;
}
