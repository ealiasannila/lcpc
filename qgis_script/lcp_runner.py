from qgis.core import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from subprocess import call
import os.path
import inspect


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

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))


call([path+"/lcpc", cost_surface, target_points,start_point, friction_field, "-o", output_path, "-p", output_points, "-d", str(max_distance), "--driver", output_driver, "-a", algorithm, "--overwrite"])
