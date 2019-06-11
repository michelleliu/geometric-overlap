#!/usr/bin/env python3

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import sys
import math
import sympy
from ase import geometry
import time

DIM = 3

def distance(x1, y1, z1, x2, y2, z2):
    dx = 1.0*(x1-x2)
    dy = 1.0*(y1-y2)
    dz = 1.0*(z1-z2)
    dist = math.sqrt(dx*dx+dy*dy+dz*dz)
    return dist

def distance_pbc(a1, a2, cell, pbc):
    D, D_len = geometry.get_distances(a1, a2, cell=cell, pbc=pbc)
    cell_trans = [0, 0, 0]
    return [D[0,0,0],D[1,1,1],D[2,2,2]], D_len[0, 0], cell_trans
    #return D, D_len[0, 0], cell_trans

def shift_point(point, delta, cell):
    shifted_point = point + np.dot(delta, cell)
    return shifted_point

def cartn_to_frac(point,invcell):
    frac = np.dot(np.transpose(invcell),point)
    return frac

def frac_to_cartn(point,cell):
    cart = np.dot(np.transpose(cell),point)
    return cart

def get_minimum_image(r0,mid):
    dr = np.zeros(DIM)
    retR = np.zeros(DIM)
    for i in range(DIM):
        dr[i] = r0[i] - mid[i]
        retR[i] = mid[i]+dr[i]
    return retR, dr

# how to run
helpline="to run:\
        \n    projection_pld_slice_PBC.py FOLDER FILENAME AXIS CUTRAD RADFILE x1 y1 z1 x2 y2 z2\
        \n    projection_pld_slice_PBC.py FOLDER FILENAME AXIS CUTRAD RADFILE x1 y1 z1 x2 y2 z2 a b c dp1 dp2 dp3 alpha beta gamma\
        \n    e.g.:\
        \n    ./projection_pld_slice_PBC.py . HKUST-1.xyz a 10 radii.txt 2.66837 13.17150 13.17150  9.85575 13.17150 13.17150\
        \n    ./projection_pld_slice_PBC.py . HKUST-1.xyz a 10 radii.txt 2.66837 13.17150 13.17150  9.85575 13.17150 13.17150 0 0 0 26.343 26.343 26.343 90 90 90\
        "
errorline="incorrect argument"
if (len(sys.argv))==2:
    if (sys.argv[1]=='-h'):
        print(helpline)
        exit(0)
    else:
        print(errorline)
        exit(0)

# read in arguments
folder = sys.argv[1]
structure = sys.argv[2]
axis = sys.argv[3]
cutrad = float(sys.argv[4])
cutrad1 = cutrad*math.sqrt(2)
cutdist_plane = 5.3
radfile = sys.argv[5] # radius file
x1 = float(sys.argv[6])
y1 = float(sys.argv[7])
z1 = float(sys.argv[8])
x2 = float(sys.argv[9])
y2 = float(sys.argv[10])
z2 = float(sys.argv[11])

# read in arguments
deltapos = np.zeros(DIM,dtype=int) # change in unit cell vector
allH = False
if (len(sys.argv) >= 21):
    deltapos[0] = float(sys.argv[12])
    deltapos[1] = float(sys.argv[13])
    deltapos[2] = float(sys.argv[14])
    cella = float(sys.argv[15])
    cellb = float(sys.argv[16])
    cellc = float(sys.argv[17])
    alpha = float(sys.argv[18])
    beta  = float(sys.argv[19])
    gamma = float(sys.argv[20])
if (len(sys.argv) >= 22):
    allH = bool(sys.argv[21])

# set unit cell
cell = geometry.cellpar_to_cell([cella, cellb, cellc, alpha, beta, gamma])
invcell = inv(cell)
print(cell)
print(invcell)

# Set centers for wrap_positions
# for when you want to specify how to wrap atoms which are
# positioned exactly on the unit cell edge.
# This is useful when, for example, node2 is shifted (1,0,0)
# and node1 has x position that lies right on UC boundary,
# so you want it to wrap to the higher value (e.g., BAZJET).
centercorrection=np.identity(3)*0.00001
center1 = np.full(DIM,0.5) + np.dot(centercorrection,deltapos)
#center2 = np.full(DIM,0.5) - np.dot(centercorrection,deltapos)

# Wrap node positions so they start in the same unit cell
# because zeo++ sometimes gives position of node1 in the
# extended unit cell.
# (Not sure this is also happening with node2, but wrap just in case)
#node1 = [x1,y1,z1]
# wrap node 1
node1 = [x1,y1,z1]
node2 = [x2,y2,z2]
print("node1: ",node1)
print("node2: ",node2)
node1 = geometry.wrap_positions([node1], cell, center=center1)[0]
node1_frac = cartn_to_frac(node1, invcell)
# wrap node 2 to same image as node 1
node2 = geometry.wrap_positions([node2], cell, center=node1_frac)[0]

# translate node2 if across unit cell
#node2 = shift_point(node2, deltapos, cell)
print("wrapped node1: ",node1)
print("wrapped node2: ",node2)

# compute midpoint and vectors
mid = np.zeros(DIM) # FIXME declare array of length DIM
edgevec = np.zeros(DIM)
edgelength = 0.0
normvec = np.zeros(DIM)
for i in range(DIM):
    mid[i] = 0.5*(node1[i]+node2[i])
    edgevec[i] = node2[i]-node1[i]
    edgelength += edgevec[i]*edgevec[i]
edgelength = math.sqrt(edgelength) # length of edgevec
for i in range(DIM):
    normvec[i] = edgevec[i]/edgelength
wrappedmidpt = geometry.wrap_positions([mid], cell)[0]
print("wrapped midpoint: ",wrappedmidpt)
print("cut: ",cutrad,"midpoint: ",mid[0],mid[1],mid[2],"vec: ",edgevec[0],edgevec[1],edgevec[2],"normvec: ",normvec[0],normvec[1],normvec[2])

# declare sympy Point3D from UNWRAPPED midpoint
# FIXME: not sure if this is the best thing to do though? should not matter if it is wrapped
midpoint = sympy.Point3D(mid[0],mid[1],mid[2])
nodepoint1 = sympy.Point3D(node1[0],node1[1],node1[2])
nodepoint2 = sympy.Point3D(node2[0],node2[1],node2[2])
vecline = sympy.geometry.Line(nodepoint1,nodepoint2)

# declare sympy Plane
plane = sympy.Plane(sympy.Point3D(mid[0],mid[1],mid[2]),normal_vector=normvec)

# declare radii dictionary
radii = dict()

# declare containers to store atoms within cutoff
atomname = [] # atom names
point3D = [] # atom centers
projection3D = [] # projection of atom centers onto 3D plane
projection2D = [] # parameterization of atom centers onto a 2D plane
radius = [] # radius of atoms
indices = []
allatomnames = []
allatompos = []
allatomnames_supercell = []
allatompos_supercell = []

# read radius file into dictionary
r = open(radfile)
for line in r:
    s = line.split()
    radii[s[0]] = float(s[1])

# open structure file
filename=folder+"/"+structure+".xyz"
f = open(filename)

# read structure file
numatoms = int(f.readline())
f.readline()

index=0
mid_frac = cartn_to_frac(mid,invcell)
supercell=2*cell
mid_frac_supercell = cartn_to_frac(mid, inv(supercell))

print("reading in atoms and projecting onto 3D plane...")
for i in range(numatoms):
    line = f.readline()
    s = line.split()
    atomtype = s[0]
    x = float(s[1])
    y = float(s[2])
    z = float(s[3])
    allatomnames.append(atomtype) # FIXME make supercell
    allatompos.append([x,y,z])

for i in range(numatoms):
    pos = allatompos[i]
    atomtype = allatomnames[i]
    for repA in range(2):
        for repB in range(2):
            for repC in range(2):
                if(repA+repB+repC > 0):
                    rep = [repA, repB, repC]
                    reppos = shift_point(pos, rep, cell)
                    reppos_shifted = geometry.wrap_positions([[reppos[0],reppos[1],reppos[2]]], supercell, center=mid_frac_supercell)[0]
                    allatomnames.append(atomtype)
                    allatompos.append(reppos_shifted)


cell=supercell
invcell = inv(cell)

fdebug = open(structure+"_"+axis+"axis"+"_supercell.xyz","w")
fdebug.write(str(len(allatomnames)))
fdebug.write("\n")
for index in range(len(allatomnames)):
    pos = allatompos[index]
    atomtype = allatomnames[index]
    fdebug.write("\n")
    fdebug.write(atomtype+" "+str(pos[0])+" "+str(pos[1])+" "+str(pos[2]))
fdebug.close()

numsupercellatoms=len(allatomnames)
for index in range(numsupercellatoms):
    if (index % round(numsupercellatoms/10) == 0):
        print(str(index)+" of "+str(numsupercellatoms))
    pos = allatompos[index]
    atomtype = allatomnames[index]
    r0_frac = cartn_to_frac(pos, invcell)
    rij_frac, drij_frac = get_minimum_image(r0_frac, mid_frac)
    drij = frac_to_cartn(drij_frac, cell)
    rij = frac_to_cartn(rij_frac, cell)
    rij_sym = sympy.Point3D(rij[0], rij[1], rij[2])

    plane_dist = plane.distance(rij_sym)
    #pbcdist = math.sqrt(drij[0]*drij[0]+drij[1]*drij[1]+drij[2]*drij[2])
    #if pbcdist < cutrad1:
    if plane_dist < cutdist_plane: # PBC search for closest atoms in all images, write coords of closest image
        cyl_dist = vecline.distance(rij_sym)
        if cyl_dist < cutrad:
            r0_sym = sympy.Point3D(pos[0],pos[1],pos[2])
            p = plane.projection(rij_sym)
            atomname.append(atomtype)
            point3D.append(rij_sym)
            projection3D.append(p)
            if(allH):
                radius.append(radii["H"])
            else:
                radius.append(radii[atomtype])
            indices.append(index)

f.close()
numpoints = len(atomname)

# project points onto 2D plane
p0 = projection3D[0]
a1 = [p0.x-mid[0], p0.y-mid[1], p0.z-mid[2]]
a2 = np.cross(a1, normvec)# in plane and perpendicular to a1
# normalize vectors
norm=math.sqrt(a1[0]*a1[0]+a1[1]*a1[1]+a1[2]*a1[2])
a1_norm = [a1[0]/norm, a1[1]/norm, a1[2]/norm]
norm=math.sqrt(a2[0]*a2[0]+a2[1]*a2[1]+a2[2]*a2[2])
a2_norm = [a2[0]/norm, a2[1]/norm, a2[2]/norm]

foutxyz=open(structure+"_"+axis+"axis"+"_cut"+str(cutrad)+".xyz","w")
#foutproj=open(structure+"_"+axis+"axis"+"_cut"+str(cutrad)+"_proj.xyz","w")

print("projecting points onto 3D plane...")
foutxyz.write(str(numpoints))
#foutproj.write(str(numpoints))
foutxyz.write("\n")
#foutproj.write("\n")
foutxyz.write("name xw yw zw xp yp zp x2 y2 r id id_cif\n")
#foutproj.write("name xp yp zp\n")
for i in range(numpoints):
    p = projection3D[i] - midpoint
    param1 = np.dot(p, a1_norm)
    param2 = np.dot(p, a2_norm)
    #print(atomname[i], float(point3D[i][0]), float(point3D[i][1]), float(point3D[i][2]), float(projection3D[i].x), float(projection3D[i].y), float(projection3D[i].z), param1, param2, radius[i], i, indices[i])
    foutxyz.write(str(atomname[i])+" "+str(float(point3D[i][0]))+" "+str(float(point3D[i][1]))+" "+str(float(point3D[i][2]))+" "+str(float(projection3D[i].x))+" "+str(float(projection3D[i].y))+" "+str(float(projection3D[i].z))+" "+str(param1)+" "+str(param2)+" "+str(radius[i])+" "+str(i)+" "+str(indices[i]))
    #foutproj.write(str(atomname[i])+" "+str(float(projection3D[i].x))+" "+str(float(projection3D[i].y))+" "+str(float(projection3D[i].z))+" "+str(indices[i]))
    foutxyz.write("\n")
    #foutproj.write("\n")
foutxyz.close()
#foutproj.close()
