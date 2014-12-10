'''
Created on 7 Dec 2014

@author: William Shakour (billy1380)
'''
from pycsg.csg import Csg
import sys

if __name__ == '__main__':
    a = Csg.cube()
    b = Csg.sphere({'radius': 1.35})
    
    c = a.subtract(b)

    p = c.toPolygons()

    #for polygon in p:
    #    for vertex in polygon.vertices:
    #        print(vertex.pos.x, vertex.pos.y, vertex.pos.z, end=' ')

    #    print(end='_')
    vert_index = {}
    verts = []
    faces = []
    
    a = Csg.cube()
    b = Csg.sphere({'radius': 1.35})
    c = Csg.cylinder({'radius': 1})
    
    rl = sys.getrecursionlimit()
    sys.setrecursionlimit(10000)
    d = a.subtract(b).subtract(c)
    sys.setrecursionlimit(rl)

    p = d.toPolygons()
    
    for polygon in p:
        face = []
        for vertex in polygon.vertices:
            pos = (vertex.pos.x, vertex.pos.y, vertex.pos.z)
            vert_key = hash(pos)
            if vert_key not in vert_index:
                verts.append(list(pos))
                vert_index[vert_key] = len(verts) - 1
            
            face.append(vert_index[vert_key])
        faces.append(face)


    verts_out = [verts]
    faces_out = [faces]
    pass