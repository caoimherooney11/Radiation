from firedrake import *
import numpy as np
from makeCellMesh_new import makeCellMesh
from generate_keff import generate_keff
import csv
from radiation_direct import solve_direct
from calculate_error import calculateError
import time

nonlinear = True
MMS = False
calculate_error = True
generate_new_lists = True
path = "NewOutput"
domain_scale = 0.5
domain_dimensions = [domain_scale*1.0, domain_scale*1.0]
cell_scale = 0.1 
cell_mesh_name = "cell_mesh_new"
direct_mesh_name = "full_mesh"

BC = "Dirichlet"
radius = 0.25 # O(delta)
k = 1. # problems arise for larger k 
c = 1. # using c in place of lambda in thesis
xi = 0.7; L = 8.; stb = 5.67e-8
Tmin = 300 # K
def vf(u):
    return 1/(2 * pi * u)
    #return 1/(4 * pi * u**2)

delta_list = []
cond_times = []
homog_times = []
homogenisation_times = []
direct_mesh_sizes = []
direct_particle_meshes = []
direct_mesh_times = []
direct_soln_times = []
errors = []
error_times = []

warning("making cell mesh...")
t0 = time.time()
#cell_mesh = makeCellMesh(cell_mesh_name, radius, cell_scale, )
global_scale = 1.
dim = 2
periodic = True
cell_mesh = makeCellMesh(cell_mesh_name, radius, cell_scale, global_scale, dim, periodic)
t1 = time.time()
cell_mesh_time = t1-t0
warning("cell mesh made -- %f seconds" % cell_mesh_time)
cell_mesh_size = cell_mesh.num_cells()
warning("cell mesh size = %f" % cell_mesh_size)

i = 0   
for delta in [0.5, 0.25, 0.125]:#, 0.0625]:
    warning("solving for delta = %f" % delta)
    delta_list.append(delta)
    direct_scale = delta * 0.25

    def RHS(u):
        return 1e-10#10**(3/8) * delta**(-3)
    T_ = 10**(8/3) * delta**(-1/3)
    tau = Tmin/T_
    warning("tau = %f" % tau)
    
    # generate effective conductivity 
    warning("generating effective conductivity...")
    t4 = time.time()
    data_name = "Output/BigDatasets/eff_cond_%d" %i
    if nonlinear:
        if generate_new_lists:
            lists = generate_keff(data_name, cell_mesh, radius, k, tau, c, vf, nonlinear) # 1st: T, 2nd-5th: keff
        else:
            with open(data+name + ".csv", 'r') as f:
                data = csv.reader(f, delimiter=",")
                lists = [ [] for _ in range(5)]
                for row in data:
                    for i in range(5):
                        lists[i].append(float(row[i]))
    else:
        effective_conductivity = generate_keff(data_name, cell_mesh, radius, k, None, c, vf, nonlinear) 
    t5 = time.time()
    warning("effective conductivity generated")
    cond_time = t5-t4
    warning("time to generate conductivity = %f" % cond_time)
    cond_times.append(cond_time)
    
    
    # solve homogenised problem 
    t6 = time.time()
    warning("solving homogenised problem...")
    if nonlinear:
        from full_problem import solve_full
        T_homog = solve_full(domain_dimensions, lists, BC, RHS)
    else:
        from linear_full_problem import solve_full
        T_homog = solve_full(domain_dimensions, effective_conductivity, BC, RHS)
    
    t7 = time.time()
    homog_soln_time = t7-t6
    warning("homogenised problem solved -- %f seconds" % homog_soln_time)
    File(path + "homog_%d.pvd" % i).write(T_homog)
    homog_times.append(homog_soln_time)
    homogenisation_times.append(homog_soln_time + cond_time)
    
    # solve direct problem
    #import sys; sys.exit()
    direct_out = File(path + "direct_%d.pvd" %i)
    error_out = File(path + "error_%d.pvd" %i)

    (T_direct, direct_mesh_size, particle_mesh, direct_mesh_time, direct_soln_time, norm_) = solve_direct(direct_mesh_name, domain_dimensions, k, delta, radius, tau, c, vf, direct_scale, BC, RHS, nonlinear, MMS)
    direct_out.write(T_direct)
    direct_mesh_sizes.append(direct_mesh_size)
    direct_particle_meshes.append(particle_mesh)
    direct_mesh_times.append(direct_mesh_time)
    direct_soln_times.append(direct_soln_time)

    warning("direct mesh size per particle = %f" % particle_mesh)
    warning("time to mesh direct problem = %f" % direct_mesh_time)
    warning("time to solve direct problem = %f" % direct_soln_time)

    if calculate_error:
        t2 = time.time()
        (error, diff) = calculateError(T_homog, T_direct, domain_dimensions)
        t3 = time.time()
        error_time = t3-t2
        warning("error = %f" % error)
        warning("time to find error = %f" % error_time)
        error_out.write(diff)
        errors.append(error)
        error_times.append(error_time)

    i = i+1

print("** PROBLEM SOLVED **")
print("time to mesh cell problem = ", cell_mesh_time)
print("time to generate effective conductivity = ", cond_times)
print("time to solve homogenised problem = ", homog_times)
print("total homogenisation time = ", homogenisation_times)
print("direct mesh size = ", direct_mesh_sizes)
print("direct mesh size per particle = ", direct_particle_meshes)
print("time to mesh direct problem = ", direct_mesh_times)
print("time to solve direct problem = ", direct_soln_times)
print("delta = ", delta_list)

if calculate_error:
    print("time to calculate error = ", error_times)
    print("error = ", errors)
    np.savetxt(path + "/error_time.csv", error_times, delimiter = ",")
    np.savetxt(path + "/error.csv", errors, delimiter = ",")

np.savetxt(path + "/delta.csv", delta_list, delimiter=",")
np.savetxt(path + "/direct_time.csv", direct_soln_times, delimiter = ",")
np.savetxt(path + "/direct_mesh_size.csv", direct_mesh_sizes, delimiter = ",")
np.savetxt(path + "/direct_particle_mesh_size.csv", direct_particle_meshes, delimiter = ",")
np.savetxt(path + "/homog_time.csv", [homog_soln_time + cond_time], delimiter = ",")
