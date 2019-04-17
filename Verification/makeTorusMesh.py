from firedrake import *
import pygmsh
import numpy as np
from geoToMesh import geoToMesh
from torusToSquare import torusToSquare

#def getTorusMesh(eps, char, global_char, name):
test = True
char = 1.0
global_char = 0.5
name = "torus"

if test:
    geom = pygmsh.built_in.Geometry()
    #geom = pygmsh.opencascade.Geometry()
    nx = 1.0
    ext = 1.0

    circle_centres = [[-ext,0.0,0.0], [0.0,-ext,0.0], [ext,0.0,0.0], [0.0,ext,0.0]]
    circle_points = []
    
    shift = [[nx/2,0.0,0.0],[0.0,0.0,nx/2],[-nx/2,0.0,0.0],[0.0,0.0,-nx/2]]
    rotate_points = []
    rotate_points.append(circle_centres[0])
    circle_points.append(geom.add_point(np.asarray(circle_centres[0]), lcar = char))
    for i in range(len(shift)):
        char_ = char
        circle_points.append(geom.add_point(np.asarray(circle_centres[0]) + np.asarray(shift[i]), lcar = char_))
        rotate_points.append(np.asarray(circle_centres[0])+np.asarray(shift[i]))
    
    shift = [[0.0,nx/2,0.0],[0.0,0.0,nx/2],[0.0,-nx/2,0.0],[0.0,0.0,-nx/2]]
    circle_points.append(geom.add_point(np.asarray(circle_centres[1]), lcar = char))
    for i in range(len(shift)):
        circle_points.append(geom.add_point(np.asarray(circle_centres[1]) + np.asarray(shift[i]), lcar = char))
    shift = [[-nx/2,0.0,0.0],[0.0,0.0,nx/2],[nx/2,0.0,0.0],[0.0,0.0,-nx/2]]
    circle_points.append(geom.add_point(np.asarray(circle_centres[2]), lcar = char))
    for i in range(len(shift)):
        char_ = char
        circle_points.append(geom.add_point(np.asarray(circle_centres[2]) + np.asarray(shift[i]), lcar = char_))
    shift = [[0.0,-nx/2,0.0],[0.0,0.0,nx/2],[0.0,nx/2,0.0],[0.0,0.0,-nx/2]]
    circle_points.append(geom.add_point(np.asarray(circle_centres[3]), lcar = char))
    for i in range(len(shift)):
        circle_points.append(geom.add_point(np.asarray(circle_centres[3]) + np.asarray(shift[i]), lcar = char))

    circle_lines = []
    for j in range(4):
        for i in range(3):
            circle_lines.append(geom.add_circle_arc(circle_points[5*j+i+1], circle_points[5*j], circle_points[5*j+i+2]))
        circle_lines.append(geom.add_circle_arc(circle_points[5*j+4], circle_points[5*j], circle_points[5*j+1]))
    
    torus_lines = []
    for i in range(4):
        for j in range(3):
            z = circle_points[5*j+i+1].x[2]
            centre = geom.add_point(np.asarray([0.0,0.0,z]), lcar = char)
            torus_lines.append(geom.add_circle_arc(circle_points[5*j+i+1], centre, circle_points[5*(j+1)+i+1]))
        torus_lines.append(geom.add_circle_arc(circle_points[5*3+i+1], centre, circle_points[i+1]))
    
    loops = []
    for i in range(3):
        for j in range(3):
            loops.append(geom.add_line_loop([circle_lines[4*i+j], torus_lines[4*j+4+i], -circle_lines[4*i+j+4], -torus_lines[4*j+i]]))
        loops.append(geom.add_line_loop([circle_lines[4*i+3], torus_lines[i], -circle_lines[4*i+3+4], -torus_lines[4*3+i]]))
    for j in range(3):
        loops.append(geom.add_line_loop([circle_lines[4*3+j], torus_lines[4*j+7], -circle_lines[j], -torus_lines[4*j+3]]))
    loops.append(geom.add_line_loop([circle_lines[15], torus_lines[3], -circle_lines[3], -torus_lines[15]]))

    surfaces = []
    print(len(loops))
    for i in range(len(loops)):
        surfaces.append(geom.add_surface(loops[i]))
    
    # physical entities
    periodic_boundary = []
    for j in range(4):
        periodic_boundary.append(circle_lines[8+j])
    for j in range(4):
        periodic_boundary.append(torus_lines[8+j])
    periodic = geom.add_physical_line(periodic_boundary, label = "1")
    surface = geom.add_physical_surface(surfaces)
   
    code = geom.get_code()
    #name = "torus"
    #global_char = 0.1
    mesh = geoToMesh(code, name, global_char, 3)

    #return mesh
    mesh_new = torusToSquare(mesh)
    V = FunctionSpace(mesh_new, "CG", 1)
    v = Function(V)
    v.interpolate(Constant(1.0))
    File("testing.pvd").write(v)
