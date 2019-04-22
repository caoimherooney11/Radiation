from firedrake import *
import numpy as np
from makeCellMesh import makeCellMesh
from generate_keff import generate_keff
from interpolate_keff  import keff

mesh_name = "cell_mesh"
radius = 0.1
hole_scale = 0.1; global_scale = 0.5
mesh = makeCellMesh(mesh_name, radius, hole_scale, global_scale)
k = 1.

lists = generate_keff(mesh, radius, k) # 1st: T, 2nd-5th: keff



