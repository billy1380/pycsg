'''
Created on 7 Dec 2014

@author: William Shakour (billy1380)
'''
from pycsg.csg import Csg

if __name__ == '__main__':
    a = Csg.cube()
    b = Csg.sphere({'radius': 1.35})
    
    c = a.subtract(b)

    p = c.toPolygons()

    for polygon in p:
        for vertex in polygon.vertices:
            print(vertex.pos.x, vertex.pos.y, vertex.pos.z, end=' ')

        print(end='_')
    pass