import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve, generate_binary_structure
import csv
from meshpy.tet import MeshInfo, build
import meshpy.triangle as triangle
import numpy.linalg as linalg

from geomdl import fitting
from geomdl.visualization import VisMPL as vis
from scipy import optimize
from scipy.interpolate import splprep, splev

import FEM
import Interpolation
import os 
import MFS
import BIE

import time

def point_in_curve(points, poly):
    x, y = points
    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n+1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside

cwd = os.getcwd()

# Here solving Laplace's Equation with delta u = 0 and u = psi on the boundary
# Load Curve_Points
curve_pts = []
crv_x = []
crv_y = []
bnd_node = []
with open(cwd + '/curve_points_file.csv', newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        if row != [] and row != ['x_pos', 'y_pos']:
            curve_pts.append((float(row[0]), float(row[1])))
            crv_x.append(float(row[0]))
            crv_y.append(float(row[1]))
            bnd_node.append([float(row[0]), float(row[1])])

bnd_node.pop(0)
bnd_node = np.array(bnd_node)     
plt.scatter(crv_x,crv_y)
## Parametrization 
## Define the nodes as a 2D array
#nodes = np.array(curve_pts)

## Perform spline interpolation to generate a smooth curve that passes through the nodes
#tck, u = splprep(nodes.T, u=None, s=0.0, per=1)

## Evaluate the curve at a set of evenly spaced parameter values to get the x, y coordinates
#t = np.linspace(0, 1, num=100, endpoint=True)
#curve = splev(t, tck)
#print(curve)

## Plot the nodes and the curve
#plt.plot(nodes[:,0], nodes[:,1], 'ro', label='Nodes')
#plt.plot(curve[0], curve[1], 'b-', label='Curve')

def round_trip_connect(start, end):
    return [(i, i + 1) for i in range(start, end)] + [(end, start)]

def needs_refinement(vertices, area):
        bary = np.sum(np.array(vertices), axis=0) / 3
        max_area = 1 + (linalg.norm(bary, np.inf) - 1) * 0.2
        return bool(area > max_area)

t2 = time.time()
info = triangle.MeshInfo()
facets = round_trip_connect(0, len(curve_pts) - 1)
circ_start = len(curve_pts)
facets.extend(round_trip_connect(circ_start, len(curve_pts) - 1))
info.set_points(curve_pts)
info.set_facets(facets)
mesh = triangle.build(info, refinement_func=needs_refinement)

mesh_points = np.array(mesh.points)
mesh_tris = np.array(mesh.elements)
t3 = time.time()
print('time spend building the triangular mesh')
print(t3-t2)
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)
plt.show()

# let's assume f(x,y) = x^2 - y^2 on the boundary
def func(x,y):
    return x*x - y*y
t0 = time.time()
crv_x.pop(0)
crv_y.pop(0)
psi = [x*x - y*y for x, y in zip(crv_x,crv_y)] # let's assume f(x,y) = x^2 - y^2

curve_pts.pop(0)
crvpts = np.array(curve_pts)
solution = FEM.laplace_solver(mesh_points,mesh_tris, crvpts, psi)
x = []
y = []
for i in mesh_points:
    x.append(i[0])
    y.append(i[1])
 
Interpolation.interpolation(x,y,mesh_tris,solution)
t1 = time.time()
print('time spend in total with FEM')
print(t1-t0)

# Define the evaluation points
x = np.linspace(int(min(crvpts[:,0])), int(max(crvpts[:,0])), int(max(crvpts[:,0])) - int(min(crvpts[:,0])))
y = np.linspace(int(min(crvpts[:,1])), int(max(crvpts[:,1])), int(max(crvpts[:,1]) - min(crvpts[:,1])))
X, Y = np.meshgrid(x, y)
points = np.stack((X.ravel(), Y.ravel())).T
t5 = time.time()
# Solve the Laplace equation and evaluate the solution
solution_function = MFS.solve_laplace_equation(crvpts, psi)
Z = solution_function(X,Y)
t6 = time.time()
print('time spend in total with MFS')
print(t6-t5)
# Plot the solution
plt.contourf(X, Y, Z.reshape(X.shape))
plt.colorbar()

graph_nodes = []
for i in crvpts:
    graph_nodes.append(i)
graph_nodes.append(graph_nodes[0])
graph_nodes = np.array(graph_nodes)
plt.plot(graph_nodes[:, 0], graph_nodes[:, 1], 'k')
plt.show()

solution_BIE = BIE.solve_laplace_single_representation(crvpts, psi, BIE.green_func)
Z_BIE = solution_BIE(points)
plt.contourf(X, Y, Z_BIE.reshape(X.shape))
plt.colorbar()

plt.plot(graph_nodes[:, 0], graph_nodes[:, 1], 'k')
plt.show()

TrueZ = func(X,Y)
plt.contourf(X, Y, TrueZ.reshape(X.shape))
plt.colorbar()

plt.plot(graph_nodes[:, 0], graph_nodes[:, 1], 'k')
plt.show()

plt.imshow(abs(TrueZ.reshape(X.shape)- Z.reshape(X.shape)))
plt.plot(graph_nodes[:, 0] - int(min(crvpts[:,0])), graph_nodes[:, 1] - int(min(crvpts[:,1])), 'k')
plt.title("MFS Error Analysis")
plt.colorbar()
plt.show()

plt.imshow(abs(TrueZ.reshape(X.shape)- Z_BIE.reshape(X.shape)))
plt.plot(graph_nodes[:, 0] - int(min(crvpts[:,0])), graph_nodes[:, 1] - int(min(crvpts[:,1])), 'k')
plt.title("BIE Error Analysis")
plt.colorbar()
plt.show()

'''
solution_double_BIE = BIE.compute_double_layer(crvpts, psi, BIE.G_normal)
Z_BIE_Double = solution_double_BIE(points)
plt.contourf(X, Y, Z_BIE_Double.reshape(X.shape))
plt.colorbar()

plt.plot(graph_nodes[:, 0], graph_nodes[:, 1], 'k')
plt.show()
'''



'''
solution_BIE_double = BIE.solve_laplace_double_representation(crvpts, psi, BIE.green_func, BIE.G_normal, points)
plt.contourf(X, Y, solution_BIE_double.reshape(X.shape))
plt.colorbar()

plt.plot(graph_nodes[:, 0], graph_nodes[:, 1], 'k')
plt.show()
'''