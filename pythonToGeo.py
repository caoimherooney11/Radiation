from firedrake import * 
#from dolfin import *

def pythonToGeo(mesh_name, domain_dimensions, radius):#, number_of_holes):
    delta = 1.0
    nx = int(domain_dimensions[0]/delta)
    ny = int(domain_dimensions[1]/delta)
    #nz = int(domain_dimensions[2]/delta)
    #eps = epsilon*delta
    #L = sqrt(delta**2 - eps**2)/2
    #L = delta/2
    #r = sqrt(delta**2 + eps**2)/2
    #num = nx*ny*nz
    r = radius
    #n = number_of_holes
    dim = 3
    f = open("%s.geo" % mesh_name, "w+")

    f.write("SetFactory(\"OpenCASCADE\");\n")
    #f.write("Sphere(1) = {%f, %f, %f, %f, -Pi/2, Pi/2, 2*Pi};\n" % (delta/2, delta/2, delta/2, r))
    #f.write("Box(2) = {0, 0, 0, %f, %f, %f};\n" % (2*L, 2*L, 2*L))
    #f.write("Surface Loop(3) = {5, 6, 2, 4, 3, 7};\n")
    #f.write("BooleanIntersection{ Volume{1}; Delete; }{ Volume{2}; Delete; }\n")
    f.write("Rectangle(1) = {0.0, 0.0, 0.0, %f, %f, 0};\n" % (nx, ny))
    f.write("Circle(5) = {%f, %f, 0, %f, 0, 2*Pi};\n" % (nx/2,ny/2,r))
    f.write("Line Loop(2) = {5};\n")
    f.write("Plane Surface(2) = {2};\n")
    f.write("BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }")
    f.write("\n")
    
   # for i in range(nx):
   #     for j in range(ny):
   #         for k in range(nz):
   #             if i == 0 and j == 0 and k == 0:
   #                 print('na')
   #             else:
   #                 f.write("Translate {%f, %f, %f} {Duplicata { Volume{1}; }}\n" % (2*L*i, 2*L*j, 2*L*k))
   #                 f.write("\n")
   #     
   # f.write("BooleanFragments{")
   # for i in range(num):
   #     f.write("Volume{%i}; " % (i+1))
   # f.write("Delete;}{}\n")
   # f.write("\n")
   # 
   # f.write("Physical Volume(\"Vol\") = {")
   # for i in range(num - 1):
   #     f.write("%d," % (i+1))
   # f.write("%d};\n" % num)#

    f.write("Physical Line(\"dirichlet1\") = {9};\n")
    f.write("Physical Line(\"dirichlet2\") = {6};\n")
    f.write("Physical Line(\"neumann\") = {8, 7};\n")
    f.write("Physical Line(\"radiation\") = {5};\n")
   # 
    f.write("Physical Surface(\"Surf\") = {1};\n")

   # limit = 7*(num) - nx*ny*(nz-1) - ny*nz*(nx-1) - nx*nz*(ny-1)
   # for i in range(limit-1):
   #     f.write("%d," % (i+1))
   # f.write("%d};\n" % (limit))

    return None




