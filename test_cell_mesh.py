from firedrake import *
from makeCellMesh import makeCellMesh

mesh_name = "cell_mesh"; radius = 0.2; scale = 0.015; global_scale = 1.0
mesh = makeCellMesh(mesh_name, radius, scale, global_scale)
