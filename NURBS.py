from geomdl import BSpline
from geomdl import knotvector
# Import Matplotlib visualization module
from geomdl.visualization import VisMPL

# Create the curve instance
def createCurveInstance():
    crv = BSpline.Curve()
    return crv

# Set degree
def setDegree(crv, degree):
    crv.degree = degree

# Set control points
def setControlPts(crv, pts):
    crv.ctrlpts = pts

# Set knot vector
def setKnotVector(crv):
    crv.knotvector = knotvector.generate(crv.degree, crv.ctrlpts_size)

# Insert knot
def insertKnor(crv, val):
    crv.insert_knot(0.3)

# Set the visualization component of the curve
def vis(crv):
    crv.vis = VisMPL.VisCurve2D()

# Plot the curve
def crvPlot(crv):
    crv.render()
    
def geoSetUp(deg, pts):
    crv = BSpline.Curve()
    crv.degree = deg
    crv.ctrlpts = pts
    crv.knotvector = knotvector.generate(crv.degree, crv.ctrlpts_size)
    return crv

def crvShow(crv):
    crv.vis = VisMPL.VisCurve2D()
    crv.render()