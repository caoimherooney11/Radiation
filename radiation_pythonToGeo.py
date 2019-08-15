from firedrake import *


def pythonToGeo(mesh_name, radius, scale):
    scale2 = 1/scale
    f = open("%s.geo" % mesh_name, "w+")

    f.write("SetFactory(\"OpenCASCADE\");\n")

    f.write("Rectangle(1) = {0, 0, 0, 1, 1, 0};\n")
    f.write("Disk(2) = {0.5, 0.5, 0, %f, %f};\n" % (radius, radius))
    f.write("BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }\n")

    f.write("Physical Line(\"hole\") = {5};\n")
    f.write("Physical Line(\"outer\") = {9, 7, 6, 8};\n")
    f.write("Physical Surface(\"surf\") = {1};\n")

    f.write("Transfinite Line {9, 8, 6, 7} = %f Using Progression 1;\n" % scale2)
    f.write("Mesh.CharacteristicLengthMin=%f;\n" % scale)
    f.write("Mesh.CharacteristicLengthMax=%f;\n" % scale)

    return None
