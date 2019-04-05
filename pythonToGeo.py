from firedrake import * 
#from dolfin import *

def pythonToGeo(mesh_name, domain_dimensions, radius, holes_x, holes_y):#, holes_z):
    #delta = 0.5
    dim = len(domain_dimensions)
    r = radius
    n = holes_x * holes_y #* holes_z

    nx = int(domain_dimensions[0])
    ny = int(domain_dimensions[1])
    if dim == 3:
        nz = int(domain_dimensions[2])
        xspace = nz/holes_z
        if 2*r > zspace:
            print("Too many circles or circles too big.")
            exit()
                

    xspace = nx/holes_x
    yspace = ny/holes_y
    if (2*r >= xspace) or (2*r >= yspace):
        print("Too many circles or circles too big.")
        exit()

    f = open("%s.geo" % mesh_name, "w+")

    f.write("SetFactory(\"OpenCASCADE\");\n")
    if dim == 2:
        f.write("Rectangle(1) = {0.0, 0.0, 0.0, %f, %f, 0};\n" % (nx, ny))
        k = 0
        for i in range(holes_x):
            for j in range(holes_y):
                f.write("Circle(%d) = {%f, %f, 0, %f, 0, 2*Pi};\n" % (5+k, xspace*(1/2 + i), yspace*(1/2 + j) , r))
                f.write("Line Loop(%d) = {%d};\n" % (2+k, 5+k))
                f.write("Plane Surface(%d) = {%d};\n" % (2+k, 2+k))
                k = k + 1
        f.write("BooleanDifference{ Surface{1}; Delete; }{")
        for i in range(n):
            f.write(" Surface{%d}; " % (2+i))
        f.write(" Delete; }")
        f.write("\n")

        f.write("Physical Line(\"dirichlet1\") = {9};\n")
        f.write("Physical Line(\"dirichlet2\") = {6};\n")
        f.write("Physical Line(\"neumann\") = {8, 7};\n")
        f.write("Physical Line(\"radiation\") = {5};\n")
        f.write("Physical Surface(\"Surf\") = {1};\n")

    else:
        # need to update 3D code for multiple holes
        f.write("Box(1) = {0, 0, 0, %f, %f, %f};\n" % (delta, delta, delta))
        f.write("Sphere(2) = {%f, %f, %f, %f, -Pi/2, Pi/2, 2*Pi};\n" % (delta/2, delta/2, delta/2, r))
        f.write("Surface Loop(3) = {5, 6, 2, 4, 3, 7};\n")
        f.write("BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }\n")
        f.write("\n")

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if i == 0 and j == 0 and k == 0:
                        print('na')
                    else:
                        f.write("Translate {%f, %f, %f} {Duplicata { Volume{1}; }}\n" % (delta*i, delta*j, delta*k))
                        f.write("\n")
            
        f.write("BooleanFragments{")
        for i in range(num):
            f.write("Volume{%i}; " % (i+1))
        f.write("Delete;}{}\n")
        f.write("\n")

        f.write("Physical Volume(\"Vol\") = {")
        for i in range(num - 1):
            f.write("%d," % (i+1))
        f.write("%d};\n" % num)#
    
        #f.write("Physical Surface(\"dirichlet1\") = {10};\n")
        #f.write("Physical Surface(\"dirichlet2\") = {12};\n")
        #f.write("Physical Surface(\"neumann\") = {9, 8, 11, 13};\n")
        #f.write("Physical Surface(\"radiation\") = {7};\n")
        #f.write("Physical Volume(\"Vol\") = {1};\n")

    return None




