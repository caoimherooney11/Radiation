from firedrake import *
import numpy as np
from makeCellMesh import makeCellMesh
from generate_keff import generate_keff
import csv
from full_problem import solve_full

mesh_name = "cell_mesh"
radius = 0.1
hole_scale = 0.1; global_scale = 0.5
mesh = makeCellMesh(mesh_name, radius, hole_scale, global_scale)
#k = 1.
def k(u):
    return u + 1.

# ------------- Generate effective conductivity ------------- # 
lists = generate_keff(mesh, radius, k) # 1st: T, 2nd-5th: keff
#with open("Output/BigDatasets/effective_conductivity.csv", 'r') as f:
#    data = csv.reader(f, delimiter=",")
#    lists = [ [] for _ in range(5)]
#    for row in data:
#        for i in range(5):
#            lists[i].append(float(row[i]))


# ---------------- Solve homogenised problem ---------------- #
domain_dimensions = [1.0, 1.0]
T = solve_full(domain_dimensions, lists)
