from firedrake import *
from makeCellMesh_new import makeCellMesh
from generate_keff import generate_keff
import numpy as np

i = 0
dim = 2.
mesh_name = "cell_mesh_cond"
radius = 0.25
#scale = 0.05
cell_scale = 0.1
global_scale = 0.5
k = 1.
Tmin = 300
c = 1.
def vf(u):
    return 1/((2 * u)**(dim-1) * pi)
nonlinear = True

delta = 0.5
T_ = 10**(8/3) * delta**(-1/3)
tau = Tmin/T_

i = 0
radius_list = [0.1, 0.2, 0.3, 0.4]
np.savetxt("Output/BigDatasets/Plots/radius.csv", radius_list, delimiter = ",")

for radius in radius_list:
    data_name = "Output/BigDatasets/Plots/eff_cond_%d" % i
    mesh =  makeCellMesh(mesh_name, radius, cell_scale, global_scale, dim, True)
    lists =  generate_keff(data_name, mesh, radius, k, tau, c, vf, nonlinear)
    i = i+1
