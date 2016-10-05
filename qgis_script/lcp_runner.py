from qgis.core import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from subprocess import call
from os.path import expanduser

home = expanduser("~")

#Here we define the input and outputs
#====================================
##cost_surface=vector
##friction_field=field cost_surface
##max_distance=number 0
##start_point=vector
##target_points=vector
##output_path=output file
##output_points=output file
##output_driver=string ESRI Shapefile
##algorithm=string astar


#And here is the body of the algorithm
#=======================================

progress.setInfo("See console for progress and additional information. Cancel button won't actually cancel execution of algorithm.")
call([home+"/.qgis2/processing/scripts/lcpc", cost_surface, target_points,start_point, friction_field, output_path, output_points, str(max_distance), output_driver, algorithm, "-o"])
