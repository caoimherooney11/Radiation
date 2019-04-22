from firedrake import *
import numpy as np
from scipy import interpolate

def keff(u, lists):
    matrix = np.zeros((2,2))
    l = 1
    for i in range(2):
        for j in range(2):
            temp = interpolate.interp1d(lists[0], lists[l])
            matrix[i,j] = temp(u)
            l = l + 1
    return matrix


