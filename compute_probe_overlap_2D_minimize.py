#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import shapely as sh
import shapely.geometry.point as shp
from matplotlib.patches import Polygon
from os.path import basename
from scipy.optimize import minimize
from scipy.optimize import Bounds

def create_ellipse(center, lengths, angle=0):
    """
    create a shapely ellipse. adapted from
    https://gis.stackexchange.com/a/243462
    """
    circ = shp.Point(center).buffer(1)
    ell = sh.affinity.scale(circ, (lengths[0]), (lengths[1]))
    ellr = sh.affinity.rotate(ell, angle)
    return ellr

def plot_ellipse(x, y, l1, l2, angle, ax):
    ellipse = create_ellipse((x, y), (l1, l1), angle)
    verts = np.array(ellipse.exterior.coords.xy)
    patch = Polygon(verts.T,color = 'red', alpha = 0.5)
    return patch

def plot_circle(x, y, r, ax):
    ellipse = create_ellipse((x, y), (r, r), 0)
    verts = np.array(ellipse.exterior.coords.xy)
    patch = Polygon(verts.T,color = 'red', alpha = 0.5)
    ax.add_patch(patch)
    return patch

def compute_circle_intersect(x, y, r, ellipse_probe):
    circle = create_ellipse((x, y), (r, r), 0)
    intersect = circle.intersection(ellipse_probe)
    return intersect.area

def plot_circle_intersect(x, y, r, ellipse_probe, ax):
    circle = create_ellipse((x, y), (r, r), 0)
    verts = np.array(circle.exterior.coords.xy)
    patch = Polygon(verts.T,color = 'red', alpha = 0.5)
    ax.add_patch(patch)
    intersect = circle.intersection(ellipse_probe)
    return intersect.area

def compute_total_intersect(variables, r1, r2, x, y, r, num_atoms):
    xshift=variables[0]
    yshift=variables[1]
    angle=variables[2]
    intersect = 0.0
    ellipse1 = create_ellipse((xshift,yshift), (r1,r2), angle)
    for i in range(num_atoms):
        area = compute_circle_intersect(x[i], y[i], r[i], ellipse1)
        intersect += area
    return intersect

# read argument
xshift = 0.0
yshift = 0.0
r1=3.5
r2=5
filepath=sys.argv[1]
filename=basename(filepath)
if(len(sys.argv) >= 4):
    xshift = float(sys.argv[2])
    yshift = float(sys.argv[3])
if(len(sys.argv) >= 6):
    r1 = float(sys.argv[4])
    r2 = float(sys.argv[5])

fig,ax = plt.subplots()

## these next few lines are pretty important because
## otherwise your ellipses might only be displayed partly
## or may be distorted
ax.set_xlim([-20,20])
ax.set_ylim([-20,20])
ax.set_aspect('equal')
ax.set_xlabel("distance [A]")
ax.set_ylabel("distance [A]")

minIntersect=100
angleMin = -1

# read in framework atoms
x = []
y = []
r = []
f=open(filepath)
num_atoms=int(f.readline())
f.readline()
for i in range(num_atoms):
    line=f.readline()
    s=line.split()
    x.append(float(s[7]))
    y.append(float(s[8]))
    r.append(float(s[9]))
    plot_circle(x[i], y[i], r[i], ax)
f.close()

angleMin = 0
xsMin = 0
ysMin = 0

for angle in np.linspace(0.0,180.0,100):
#    for xshift in np.linspace(-0.5,0.5,10):
#        for yshift in np.linspace(-0.5,0.5,10):
    intersect = compute_total_intersect([xshift, yshift, angle], r1, r2, x, y, r, num_atoms)
    if (intersect < minIntersect):
        minIntersect = intersect
        angleMin = angle
        xsMin = xshift
        yxMin = yshift
bounds=Bounds([-1, -1, 0],[1, 1, 180])
v0=np.array([xsMin,ysMin,angleMin])
res = minimize(compute_total_intersect, v0, args=(r1, r2, x, y, r, num_atoms),
        method='L-BFGS-B', bounds=bounds, options={'ftol': 1e-2, 'disp': False})
xsMin = res.x[0]
ysMin = res.x[1]
angleMin = res.x[2]
minIntersect = compute_total_intersect(res.x, r1, r2, x, y, r, num_atoms)

## plot probe ellipse in blue
ellipse1 = create_ellipse((xsMin,ysMin), (r1,r2), angleMin)
verts1 = np.array(ellipse1.exterior.coords.xy)
patch1 = Polygon(verts1.T, color = 'blue', alpha = 0.5)
ax.add_patch(patch1)

print(filename[:-4]+' '+str(minIntersect)+' '+str(angleMin)+' '+str(xsMin)+' '+str(ysMin))

plt.savefig(filename[:-4]+"_minintersect.png")
