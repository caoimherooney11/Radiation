from firedrake import *
import copy
import numpy as np


def calculateError(f_homog, f_direct, domain_dimensions):#, name):
    warning("calculating error...")
    # interpolate f_homog into f_direct.function_space()
    dofs_direct = Function(FunctionSpace(f_direct.ufl_domain(), VectorElement(f_direct.function_space().ufl_element())))
    dofs_direct.interpolate(SpatialCoordinate(f_direct.ufl_domain()))
    target = f_direct.copy(deepcopy=True)
    targetvec = target.vector().get_local()
    dofsvec = dofs_direct.vector()[:, :].copy()
    for i in range(len(dofsvec[:,:])):
        x = dofsvec[i,:] 
        x = [x[0], x[1]]
        x[0] = min(max(x[0], 0), domain_dimensions[0])
        x[1] = min(max(x[1], 0), domain_dimensions[1])
        targetvec[i] = f_homog.at(x, tolerance=1e-6)
    temp = target.vector()
    temp.set_local(targetvec)
    error = f_direct.copy(deepcopy=True)
    error -= target
    warning("error calculated = %f" % norm(error))
    return norm(error), error

