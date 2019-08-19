from firedrake import *
import numpy as np
from scipy import interpolate

def keff(u, lists, entry=False, derivative=False):
    
    bbox = [-1e-6, 1+1e-6]
    if entry:
        temp = interpolate.InterpolatedUnivariateSpline(lists[0], lists[int(entry)], bbox=bbox)
        if derivative:
            out = temp.derivatives(u)[int(derivative)]
        else:
            out = temp(u)

    else:
        matrix = np.zeros((2,2))
        l = 1
        for i in range(2):
            for j in range(2):
                #temp = interpolate.interp1d(lists[0], lists[l])
                temp = interpolate.InterpolatedUnivariateSpline(lists[0], lists[l], bbox=bbox)
                if derivative:
                    matrix[i,j] = temp.derivatives(u)[int(derivative)] 
                else:
                    matrix[i,j] = temp(u)
                l = l + 1
        out = matrix

    return out
    


