from firedrake import * 
#from dolfin import *

def pythonToGeo(mesh_name, radius, scale):
    f = open("%s.geo" % mesh_name, "w+")

    f.write("SetFactory(\"OpenCASCADE\");\n")
    f.write("Torus(1) = {0.0, 0.0, 0.0, 1.0, 0.5, 2*Pi};\n")
    f.write("Sphere(2) = {-0.5, 0.0, 0.0, %f, -Pi/2, Pi/2, 2*Pi};\n" % (radius * 3))
    f.write("Delete{Volume{2}; Volume{1};}\n")
    f.write("BooleanDifference{Surface{1}; Delete;}{Surface{2}; Delete;}\n")
    f.write("Delete { Surface{2}; }\n")
    f.write("\n")

    f.write("Physical Line(\"hole\") = {3};\n")
    f.write("Physical Line(\"line1\") = {1};\n")
    f.write("Physical Line(\"line2\") = {2};\n")
    f.write("Physical Surface(\"surf\") = {1};\n")
    f.write("\n")

    f.write("Characteristic Length{2} = %f;\n" % scale)
    return None




