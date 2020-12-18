#%% matplotlib inline
import time
from matplotlib.pyplot import subplots
from pyFAI.gui import jupyter
import pyFAI
import fabio
from pyFAI.test.utilstest import UtilsTest
from pyFAI.calibrant import CALIBRANT_FACTORY
from pyFAI.goniometer import SingleGeometry
print(pyFAI.version)
start_time = time.perf_counter()
# %%
# In this example, we will re-use one of the image used int the test-suite
filename = '/Users/rbelisle/Desktop/LaB6_det315_3s_01221009_0001.tif'
frame = fabio.open(filename).data
# and now display the image
ax = jupyter.display(frame)
# %%
# This allow to measure approximatively the position of the beam center ...
x = 1500 # x-coordinate of the beam-center in pixels
y = 2750 # y-coordinate of the beam-center in pixels
d =  315 # This is the distance in mm (unit used by Fit2d)
wl = 0.9762e-10 # The wavelength is 1 Å
# %%
# Definition of the detector and of the calibrant:
pilatus = pyFAI.detectors.Detector(pixel1 = 73.242e-6, pixel2 = 73.242e-6, max_shape=(.225, .225))
behenate = CALIBRANT_FACTORY("LaB6")
behenate.wavelength = wl
behenate
# %%
# Set the guessed geometry
initial = pyFAI.geometry.Geometry(detector=pilatus, wavelength=wl)
initial.setFit2D(d,x,y)
initial
# %% The SingleGeometry object (from goniometer) allows to extract automatically ring and calibrate
sg = SingleGeometry("demo", frame, calibrant=behenate, detector=pilatus, geometry=initial)
sg.extract_cp(max_rings=5)
# %% Control point and rings do not overlap well initially (this was a guessed geometry)
jupyter.display(sg=sg)
# %% Refine the geometry ... here in SAXS geometry, the rotation is fixed in orthogonal setup
sg.geometry_refinement.refine2(fix=["rot1", "rot2", "rot3", "wavelength"])
sg.get_ai()

# %%
