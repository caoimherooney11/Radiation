from firedrake import *


def pythonToGeo(mesh_name, radius, scale, dim):
    scale2 = 1/scale
    f = open("%s.geo" % mesh_name, "w+")
    f.write("SetFactory(\"OpenCASCADE\");\n")

    if dim is 2:
        f.write("Rectangle(1) = {0, 0, 0, 1, 1, 0};\n")
        f.write("Disk(2) = {0.5, 0.5, 0, %f, %f};\n" % (radius, radius))
        f.write("BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }\n")

        f.write("Physical Line(\"hole\") = {5};\n")
        f.write("Physical Line(\"outer\") = {9, 7, 6, 8};\n")
        f.write("Physical Surface(\"surf\") = {1};\n")

        f.write("Transfinite Line {9, 8, 6, 7} = %f Using Progression 1;\n" % scale2)
        f.write("Mesh.CharacteristicLengthMin=%f;\n" % scale)
        f.write("Mesh.CharacteristicLengthMax=%f;\n" % scale)

    elif dim is 3:
        f.write("Box(1) = {0, 0, 0, 1, 1, 1};\n")
        f.write("Sphere(2) = {0.5, 0.5, 0.5, %f, -Pi/2, Pi/2, 2*Pi};\n" % radius)
        f.write("Surface Loop(3) = {5, 6, 2, 4, 3, 7};\n")
        f.write("BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }\n")
        f.write("\n")

        f.write("Physical Surface(\"hole\") = {7};\n")
        f.write("Physical Surface(\"outer\") = {8, 9, 10, 11, 12, 13};\n")
        f.write("Physical Volume(\"vol\") = {1};\n")

        f.write("Transfinite Surface {8, 9, 10, 11, 12, 13};")
        scale_min = scale
        scale_max = scale#*20   
        f.write("Mesh.CharacteristicLengthMin=%f;\n" % scale_min)
        f.write("Mesh.CharacteristicLengthMax=%f;\n" % scale_max)

        f.write("Field[1] = Ball;\n")
        radius_ball = radius*1.5
        f.write("Field[1].Radius = %f;\n" % radius_ball)
        f.write("Field[1].VIn = %f;\n" % scale_min)
        f.write("Field[1].VOut = %f;\n" % scale_max)
        f.write("Field[1].XCenter = 0.5;\n")
        f.write("Field[1].YCenter = 0.5;\n")
        f.write("Field[1].ZCenter = 0.5;\n")
        f.write("Background Field = 1;\n")

    return None
