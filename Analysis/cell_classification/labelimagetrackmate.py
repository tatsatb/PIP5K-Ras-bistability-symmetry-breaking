# This is a Jython(Python) script that uses the TrackMate to analyze and track cells in a hyperstack image in an
# automated, headless way. Since this is a Jython script, it is not meant to be run in a standard Python environment,
# but on a Java virtual machine (JVM) with Jython installed.

# This script can be executed in Fiji/ImageJ environment. For details on how to run Jython programs in Fiji/ImageJ,
# please refer to the official documentation: https://imagej.net/scripting/jython/

# This program was craeted and tested on Fiji/ImageJ 1.54p with TrackMate v7.14.0 installed (OS: Ubuntu 22.04 LTS).



import sys

from ij import IJ
from ij import WindowManager

from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.detection import MaskDetectorFactory  
from fiji.plugin.trackmate.tracking.kdtree import NearestNeighborTrackerFactory 
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.gui.displaysettings.DisplaySettings import TrackMateObject
from fiji.plugin.trackmate.features.track import TrackIndexAnalyzer

import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter

from fiji.plugin.trackmate.action import LabelImgExporter
from fiji.plugin.trackmate.util import LogRecorder
from fiji.plugin.trackmate.action.LabelImgExporter.LabelIdPainting import LABEL_IS_TRACK_ID

# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')

# Get currently selected image
imp = WindowManager.getCurrentImage()

# Swap t and z if required

if imp is not None:
    nFrames = imp.getNFrames()
    nSlices = imp.getNSlices()
    
    print("Current dimensions: T=%d, Z=%d" % (nFrames, nSlices))
    
    # If t=1 and Z>1, exchange Z and T
    if nFrames == 1 and nSlices > 1:
        print("Swapping Z and T dimensions (Z=%d -> T=%d)" % (nSlices, nSlices))
        # Create a new hyperstack with exchanged dimensions
        imp = IJ.run(imp, "Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]")
        imp = WindowManager.getCurrentImage()
        print("New dimensions: T=%d, Z=%d" % (imp.getNFrames(), imp.getNSlices()))



imp.show()


#----------------------------
# Create the model object now
#----------------------------

# Some of the parameters we configure below need to have
# a reference to the model at creation. So we create an
# empty model now.

model = Model()

# Send all messages to ImageJ log window.
model.setLogger(Logger.IJ_LOGGER)



#------------------------
# Prepare settings object
#------------------------

settings = Settings(imp)

# Configure detector - We use the Strings for the keys
settings.detectorFactory = MaskDetectorFactory()
settings.detectorSettings = {
    'TARGET_CHANNEL' : 1,
    'SIMPLIFY_CONTOURS' : True,
}  

# Configure spot filters - Classical filter on quality
filter1 = FeatureFilter('QUALITY', 66, True)
settings.addSpotFilter(filter1)

# Configure tracker - We want to allow merges and fusions
settings.trackerFactory = NearestNeighborTrackerFactory()
settings.trackerSettings = settings.trackerFactory.getDefaultSettings() # almost good enough
settings.trackerSettings['LINKING_MAX_DISTANCE'] = 150.0


# Print the available settings - Python 2.7 compatible version
print("Available tracker settings:")
for key in settings.trackerSettings.keySet():
    print("  - %s: %s" % (key, settings.trackerSettings[key]))

# You can also print the full settings object
print("Full settings object:")
print(settings.trackerSettings)



# Now you can set the appropriate parameter based on what's displayed
# Instead of 'KEY_LINKING_MAX_DISTANCE', use the correct key name

# Set the maximal linking distance (in physical units)


# Add ALL the feature analyzers known to TrackMate. They will 
# yield numerical features for the results, such as speed, mean intensity etc.
settings.addAllAnalyzers()

# Configure track filters - We want to set an area filter to get rid of the
# small non-cell spots.

filter2 = FeatureFilter('Area', 1150, True)
settings.addTrackFilter(filter2)


#-------------------
# Instantiate plugin
#-------------------

trackmate = TrackMate(model, settings)

#--------
# Process
#--------

ok = trackmate.checkInput()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))

ok = trackmate.process()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))


#----------------
# Display results
#----------------

# A selection.
selectionModel = SelectionModel( model )

# Read the default display settings.
ds = DisplaySettingsIO.readUserDefault()
# Color by tracks.
ds.setTrackColorBy( TrackMateObject.TRACKS, TrackIndexAnalyzer.TRACK_INDEX )
ds.setSpotColorBy( TrackMateObject.TRACKS, TrackIndexAnalyzer.TRACK_INDEX )

displayer =  HyperStackDisplayer( model, selectionModel, imp, ds )
displayer.render()
displayer.refresh()

# Echo results with the logger we set at start:
model.getLogger().log( str( model ) )

logger = LogRecorder( Logger.VOID_LOGGER )

lblImg = LabelImgExporter()
lblImg.createLabelImagePlus(trackmate, False, False, LABEL_IS_TRACK_ID, logger).show()
