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
#include "sorted_vector.h"

#include <tr1/functional>
#include <tr1/unordered_set>

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

void readCostSurface(const char* costSurface, const char* targets, const char* start, LcpFinder* finder, const char* frictionField, const char* maxDist, OGRSpatialReference* sr) {
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

bool fileExists(const std::string& name) {
    if (FILE * file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
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

void writePoints(std::deque<const Coords*> results, std::string outputfile, const char *pszDriverName, OGRSpatialReference sr) {
    OGRSFDriver *driver;
    OGRRegisterAll();
    driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName);

    if (!testDriver(driver)) {
        driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
    }

    OGRDataSource *pointDS;

    outputfile = validateOutfile(driver, outputfile);
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

    for (const Coords* point: results) {
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

void writeOutput(std::deque<const Coords*> results, std::string outputfile, const char *pszDriverName, OGRSpatialReference sr) {

    OGRSFDriver *driver;
    OGRRegisterAll();
    driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName);

    if (!testDriver(driver)) {
        driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
    }

    OGRDataSource *pathDS;

    outputfile = validateOutfile(driver, outputfile);
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

int main(int argc, char* argv[]) {

    //"${OUTPUT_PATH}" large_sp.shp ltargets.shp lstart.shp Luokka3 output 500
    //"${OUTPUT_PATH}" testarea.shp targets.shp start.shp friction output 500
   
    
    if (argc < 2 or strcmp(argv[1], "-h") == 0) {
        std::cout << "This program is used to search for least cost paths in polygonal costsurface.\n";

        std::cout << "USAGE:\n"
                "The program requires 8 parameters to run correctly. Parameters should be in following order:\n"
                "\tname of cost surface file. (Supports formats supported by OGR)\n"
                "\tname of target points file can be any number of points. All points must be inside cost surface polygons. (Supports formats supported by OGR)\n"
                "\tname of start point file, must be single point inside one of the cost surface polygons. (Supports formats supported by OGR)\n"
                "\tname of the friction field in cost surface\n"
                "\toutput file name, if no extension is given a folder will be assumed (at least with shapefile driver)\n"
                "\tmaximum distance between nodes in cost surface. This is used to add temporary additional nodes during LCP calculation if nodes are too far apart.\n"
                "\toutput driver name (Supports drivers supported by OGR)\n"
                "\talgorithm to use in lcp search. Options are dijksra/astar. Using astar is recommended unless number of targetpoints is high compared to number of nodes.";

        exit(0);
    }

    if (argc != 10) {
        std::cout << "Invalid arguments provided see \"lcp -h\" for details\n";
        std::cout << "Correct form is\n"
                "lcpc cost_surface_file target_file start_file name_of_friction_field output_file_name max_distance_between_nodes output_drivername search_algorithm\n";
        return 0;
    }
    int alg = 1;
    if (strcmp(argv[9], "astar") == 0) {
        alg = 1;
    } else if (strcmp(argv[9], "dijkstra") == 0) {
        alg = 0;
    } else {
        std::cout << "unknown search algorithm. using A* (possible choises are \"astar\" or \"dijkstra\"\n";
    }
    OGRRegisterAll();
    LcpFinder finder{};
    std::cout << "Reading cost surface...\n";
    OGRSpatialReference sr;
    std::clock_t begin = std::clock();
    readCostSurface(argv[1], argv[2], argv[3], &finder, argv[4], argv[7], &sr);
    //readCostSurfaceDummy("testpolygon.shp", "targets.shp", "start.shp", &finder);
    double secs = double(std::clock() - begin) / CLOCKS_PER_SEC;
    std::cout << "Finished reading cost surface (took " << secs << " s). Starting LCP search...\n";
    begin = std::clock();
    std::deque<const Coords*> results = finder.leastCostPath(alg);
    secs = double(std::clock() - begin) / CLOCKS_PER_SEC;
    std::cout << "\nSearch finished (took " << secs << "s). Writing results...\n";


    writeOutput(results, argv[5], argv[8], sr);
    writePoints(results, argv[6], argv[8], sr);
    std::cout << "All done!\n";
    return 0;
}
