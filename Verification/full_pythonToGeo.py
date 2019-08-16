from firedrake import * 
#from dolfin import *

def pythonToGeo(mesh_name, domain_dimensions, delta, radius):
    #delta = 0.5
    dim = len(domain_dimensions)
    xlim = domain_dimensions[0]
    ylim = domain_dimensions[1]
    nx = int(xlim/delta)
    ny = int(ylim/delta)
    holes_x = nx
    holes_y = ny
    r = radius*delta
    n = holes_x * holes_y
    num = nx*ny

    #if dim == 3:
    #    nz = int(domain_dimensions[2])
    #    xspace = nz/holes_z
    #    if 2*r > zspace:
    #        print("Too many circles or circles too big.")
    #        exit()
                
    #xspace = nx/holes_x
    #yspace = ny/holes_y
    #if (2*r >= xspace) or (2*r >= yspace):
    #    print("Too many circles or circles too big.")
    #    exit()

    f = open("%s.geo" % mesh_name, "w+")

    f.write("SetFactory(\"OpenCASCADE\");\n")
    if dim == 2:
        f.write("Rectangle(1) = {0.0, 0.0, 0.0, %f, %f, 0};\n" % (xlim, ylim))
        #f.write("Rectangle(1) = {0.0, 0.0, 0.0, %f, %f, 0};\n" % (delta, delta))
        #f.write("Circle(5) = {%f, %f, 0, %f, 0, 2*Pi};\n" % (delta/2, delta/2, r))
        #f.write("Line Loop(6) = {5};\n")
        #f.write("Plane Surface(6) = {6};\n")
        #f.write("BooleanDifference{ Surface{1}; Delete; }{ Surface{6}; Delete; }\n")
        #f.write("\n")

        #for i in range(nx):
        #    for j in range(ny):
        #        if i == 0 and j == 0: 
        #            print('na')
        #        else:
        #            f.write("Translate {%f, %f, 0.0} {Duplicata { Surface{1}; }}\n" % (delta*i, delta*j))
        #            f.write("\n")
        #    
        #f.write("BooleanFragments{")
        #for i in range(num):
        #    f.write("Surface{%i}; " % (i+1))
        #f.write("Delete;}{}\n")
        #f.write("\n")
        #
        #f.write("Physical Surface(\"surf\") = {")
        #for i in range(num-1):
        #    f.write("%d," % (i+1))
        #f.write("%d};\n" % (num))

        k = 0
        for i in range(nx): 
            for j in range(ny):
                f.write("Circle(%d) = {%f, %f, 0, %f, 0, 2*Pi};\n" % (5+k, (2*i+1)*delta/2, (2*j+1)*delta/2, r))
                f.write("Line Loop(%d) = {%d};\n" % (2+k, 5+k))
                f.write("Plane Surface(%d) = {%d};\n" % (2+k, 2+k))
                k = k + 1

        f.write("BooleanDifference{ Surface{1}; Delete; }{")
        for i in range(n):
            f.write(" Surface{%d}; " % (2+i))
        f.write(" Delete; }")
        f.write("\n")

        for i in range(n):
            f.write("Physical Line(\"radiation_%d\") = {%d};\n" % (i+1, i+5))
        f.write("Physical Line(\"dirichlet1\") = {%d};\n" % (n+8))
        f.write("Physical Line(\"dirichlet2\") = {%d};\n" % (n+5))
        f.write("Physical Line(\"neumann\") = {%d, %d};\n" % (n+6, n+7))
        f.write("Physical Surface(\"Surf\") = {1};\n")

    #else:
    #    # need to update 3D code for multiple holes
    #    f.write("Box(1) = {0, 0, 0, %f, %f, %f};\n" % (delta, delta, delta))
    #    f.write("Sphere(2) = {%f, %f, %f, %f, -Pi/2, Pi/2, 2*Pi};\n" % (delta/2, delta/2, delta/2, r))
    #    f.write("Surface Loop(3) = {5, 6, 2, 4, 3, 7};\n")
    #    f.write("BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }\n")
    #    f.write("\n")

    #    for i in range(nx):
    #        for j in range(ny):
    #            for k in range(nz):
    #                if i == 0 and j == 0 and k == 0:
    #                    print('na')
    #                else:
    #                    f.write("Translate {%f, %f, %f} {Duplicata { Volume{1}; }}\n" % (delta*i, delta*j, delta*k))
    #                    f.write("\n")
    #        
    #    f.write("BooleanFragments{")
    #    for i in range(num):
    #        f.write("Volume{%i}; " % (i+1))
    #    f.write("Delete;}{}\n")
    #    f.write("\n")

    #    f.write("Physical Volume(\"Vol\") = {")
    #    for i in range(num - 1):
    #        f.write("%d," % (i+1))
    #    f.write("%d};\n" % num)#
    
        #f.write("Physical Surface(\"dirichlet1\") = {10};\n")
        #f.write("Physical Surface(\"dirichlet2\") = {12};\n")
        #f.write("Physical Surface(\"neumann\") = {9, 8, 11, 13};\n")
        #f.write("Physical Surface(\"radiation\") = {7};\n")
        #f.write("Physical Volume(\"Vol\") = {1};\n")

    return None




