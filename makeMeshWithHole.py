from firedrake import *
from makeCellMesh import makeCellMesh

mesh_name = "cell_mesh"
hole_scale = 0.1
global_scale = 0.5
radius = 0.2
mesh = makeCellMesh(mesh_name, radius, hole_scale, global_scale)
