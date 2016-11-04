#!/usr/bin/env python

# simple example for pyGRASS usage: raster processing via Modules approach
# Read GeoTIFF directly, write out GeoTIFF directly
import os
import tempfile
import grass.script as grass
import numpy as np
from sets import Set
from subprocess import call
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.raster import RasterRow

from grass.pygrass.vector import VectorTopo
from grass.pygrass.vector.geometry import Point


path = '/home/elias/Documents/lcp_test_data/ds5/'
cs = 'cs.shp'
vpath = 'visi_path.shp'
vpoints = 'v_points.shp'
vpmap = 'v_points'
start_shp = 'start.shp'
target_shp = 'targets.shp'
csmap = 'cs'
crmap = 'cr'
cumul16 = 'cumul16'
cumul8 = 'cumul8'
dir8 = 'dir8'
dir16 = 'dir16'
pr8 = 'path_ras8'
pr16 = 'path_ras16'
pv8 = 'path_vec8'
pvc8 = 'path_vec_cost8'
pv16 = 'path_vec16'
pvc16 = 'path_vec_cost16'
nogo16 = 'path_nogo16'
nogo8 = 'path_nogo8'
startmap = 'start'
targetmap = 'targets'
trg_c8 = 'cost_8'
trg_c16 = 'cost_16'
trg_cv = 'cost_visi'
fcol = 'friction'
mapset = 'ds5'

result_file = 'results.csv'

def importCs():
    v.in_ogr(path + cs, output=csmap, overwrite = True)
def rasterizeCs(reso):
    g.region(res=reso)
    v.to_rast(csmap, output = 'cr', use='attr',attribute_column=fcol, overwrite = True)
def getRandomGridCoords(n, accessible):
    c = grass.region()
    rows = c['rows']
    cols = c['cols']
    nc = c['n']
    wc = c['w']
    ns = c['nsres']
    ew = c['ewres']
    
    if(rows*cols<n):
        n = rows*cols
    rand_rows = np.random.randint(0,rows, n)
    rand_cols = np.random.randint(0,cols, n)
    if accessible:
        rvals = RasterRow(crmap, mapset)
        rvals.open('r')
        
        for i in xrange(0,n):
            for j in xrange(0,n):
                while rvals[rand_rows[i]][rand_cols[j]] < 0 or rvals[rand_rows[i]][rand_cols[j]]>=999:
                    rand_rows[i] = np.random.randint(0,rows)
                    rand_cols[i] = np.random.randint(0,cols)
        rvals.close()
    
    return [Point(wc + rand_cols[i]*ew + ew/2, nc - rand_rows[i]*ns - ns/2) for i in xrange(0,n)]
def createRandomPoints(n,mapname, accessible):
    cols = [(u'cat','INTEGER PRIMARY KEY'),(trg_cv, 'DOUBLE'),(trg_c8, 'DOUBLE'),(trg_c16, 'DOUBLE')]
    points = getRandomGridCoords(n, accessible)
    with VectorTopo(mapname, mode='w', tab_cols=cols, overwrite=True) as vmap:
        # save the point and the attribute
        for i in xrange(0,n):
            vmap.write(points[i],(-2,))
        # save the changes to the database
        vmap.table.conn.commit()

def lcp_rast(knights):
    if(knights):
        r.cost(crmap, output=cumul16, outdir=dir16, start_points=startmap, stop_points=targetmap, overwrite = True,flags = 'k')
        v.what_rast(targetmap,raster=cumul16, column=trg_c16)             
        r.drain( cumul16,direction=dir16, output=pr16, drain=pv16, start_points=targetmap, flags = 'dc', overwrite = True)
        v.overlay(ainput=pv16, binput=csmap, operator='and', output=pvc16,overwrite=True) 
        v.overlay(ainput=pv16, binput=csmap, operator='not', output=nogo16,overwrite=True)
    else:
        r.cost(crmap, output=cumul8, outdir=dir8, start_points=startmap, stop_points=targetmap, overwrite = True)
        v.what_rast(targetmap,raster=cumul8, column=trg_c8)
        r.drain( cumul8,direction=dir8, output=pr8, drain=pv8, start_points=targetmap, flags = 'dc', overwrite = True)
        v.overlay(ainput=pv8, binput=csmap, operator='and', output=pvc8,overwrite=True)
        v.overlay(ainput=pv8, binput=csmap, operator='and', output=nogo8,overwrite=True)

def exportPoints():
    v.out_ogr(targetmap, output=path + start_shp, format='ESRI_Shapefile', overwrite=True)
    v.out_ogr(startmap, output=path + targets_shp, format='ESRI_Shapefile', overwrite=True)

def lcp_visi():
    max_distance = 1000
    call(["lcpc", path + cs, path+target_shp, path+start_shp, fcol, "-o", path + vpath, "-p", path + vpoints, "-d", str(max_distance), "--driver", 'ESRI Shapefile', "-a", 'astar', "--overwrite"])
    
def resultsToTargetPoints():
    v.in_ogr(path + vpoints, output=vpmap, overwrite = True)
    v.distance(targetmap, to=vpmap, dmax=0, upload='to_attr', column=trg_cv, to_column='cost')

def sumRasterPath(pvc,nogo):
    ng_path = VectorTopo(nogo)
    ng_path.open('r')
    invalid = Set([])
    for seg in ng_path:
        invalid.add(seg.attrs['a_cat'])
    ng_path.close()
    path = VectorTopo(pvc)
    path.open('r')
    actual_costs = {}
    for seg in path:
        if seg.attrs['a_cat'] in invalid:
            continue
        if actual_costs.has_key(seg.attrs['a_cat']):
            actual_costs[seg.attrs['a_cat']]+= seg.attrs['b_friction']*seg.length()
        else:
            actual_costs[seg.attrs['a_cat']] = seg.attrs['b_friction']*seg.length()
    path.close()
    return actual_costs

def exportResultTable(cell_size):
    output = [['Cat', 'visi_cost', '8_prop', '16_prop', '8_act', '16_act']]
    raster_costs8 = sumRasterPath(pvc8, nogo8)
    raster_costs16 = sumRasterPath(pvc16, nogo16)
    trg = VectorTopo(targetmap)
    trg.open('r')
    for t in trg:
        outr = [str(t.attrs['cat'])]
        print 'cat:', t.attrs['cat']
        print 'visi cost:', t.attrs['cost']
        if t.attrs['cost']>=0:
            outr.append(str(t.attrs['cost']))
        else:
            outr.append('nf')
        if t.attrs[trg_c8]:
            print '8 proposed:',t.attrs[trg_c8]*cell_size
            outr.append(str(t.attrs[trg_c8]*cell_size))
        else:
            outr.append('nf')
            print 'no 8 path found'
        if t.attrs[trg_c16]:
            print '16 proposed:',t.attrs[trg_c16]*cell_size
            outr.append(str(t.attrs[trg_c16]*cell_size))
        else:
            print 'no 16 path found'
            outr.append('nf')            
        if raster_costs8.has_key(t.attrs['cat']):
            print '8 actual:',raster_costs8[t.attrs['cat']]
            outr.append(str(raster_costs8[t.attrs['cat']]))
        else:
            print 'invalid 8 con raster path'
            outr.append('inv')
        if raster_costs16.has_key(t.attrs['cat']):
            print '16 actual:',raster_costs16[t.attrs['cat']]
            outr.append(str(raster_costs16[t.attrs['cat']]))
        else:
            print 'invalid 16 con raster path'
            outr.append('inv')            
        print ""
        output.append(outr)
    np.savetxt(path + result_file, output, delimiter = ';', fmt = '%s')
    
    

#importCs()
cellsize = 10
#rasterizeCs(cellsize)
#createRandomPoints(1, startmap, True)
#createRandomPoints(10, targetmap, False)
#lcp_rast(False)
#lcp_rast(True)
#exportPoints()
#lcp_visi()
#resultsToTargetPoints()
exportResultTable(cellsize)