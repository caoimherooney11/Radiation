from firedrake import *

nondimensional = True
nonlinear = False
mesh_name = "test_direct"
BC = "Dirichlet"
def f(u):
    return 1e-10
xi = 0.7; stb = 5.67e-8; L = 8.
delta = 0.5
radius = 0.25
Tmin = 310.15 # K  
if nondimensional:
    from radiation_direct import solve_direct
    domain_dimensions = [1.0,1.0]
    k = 1.
    T_ = (xi/(stb*delta*L))**(1/3)
    tau = Tmin/T_
    c = 1.

else:
    from dimensional_direct import solve_direct
    domain_dimensions = [L,L]
    k = xi
    c = stb
    radius = radius*L*delta
    cell_size = delta*L
    tau = cell_size
    print(cell_size)
    print(radius)


norms = []  
for i in range(3):
    scale = 2**(-i-1)
    (u_, mesh_size, particle_mesh, mesh_time, soln_time, norm_) = solve_direct(mesh_name, domain_dimensions, k, delta, radius, tau, c, scale, BC, f, nonlinear, True)
    norms.append(norm_)
    warning("error = %f" % norm_)

print(norms)
    
