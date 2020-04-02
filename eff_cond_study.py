from firedrake import *
from makeCellMesh_new import makeCellMesh
from generate_keff import generate_keff
import numpy as np
from scipy.constants import Stefan_Boltzmann    

#   Dimensional study
dimensional = True
nonlinear = True
dim = 3
mesh_name = "cell_mesh_cond"
cell_scale = 0.1
global_scale = 0.5
stb = Stefan_Boltzmann

if dimensional:
    k = 30.
    phi = 0.1 # Kris void fraction
    R = ( 3 * phi / (4 * pi) ) ** (1/3)
    tau = 0.
    radius_list = [0.01, 0.05] # Kris' radius (m)
    path = "Output/BigDatasets/Kris_comparison/"
else:
    delta = 0.5
    k = 1.
    c = 1.
    Tmin = 300
    T_ = 10**(8/3) * delta**(-1/3)
    tau = 0.#Tmin/T_
    radius_list = [0.1, 0.2, 0.3, 0.4]
    path = "Output/BigDatasets/Plots/"

np.savetxt(path + "radius.csv", radius_list, delimiter = ",")

def vf(u):
    return 1/((2 * u)**(dim-1) * pi)

i = 0

for radius in radius_list:
    if dimensional:
        deltaL = radius / R # (m)
        #L = radius / (delta * R) # (m)
        radius = R # want radius label apply to cell radius, not Kris' radius
        c = stb * deltaL / k

    data_name = path + "eff_cond_%d" % i
    warning("making mesh for R = %f" % radius)
    mesh =  makeCellMesh(mesh_name, radius, cell_scale, global_scale, dim, True)
    warning("generating conductivity for R = %f" % radius)
    lists =  generate_keff(data_name, mesh, radius, k, tau, c, vf, nonlinear, dim, dimensional)
    import IPython; IPython.embed()
    warning("anlaysis for R = %f complete" % radius)
    i = i+1
    import sys; sys.exit()
