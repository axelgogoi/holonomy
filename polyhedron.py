# Written by Axel Gogoi
# https://github.com/axelgogoi/holonomy

# Utility classes for polyhedra

import math
import networkx as nx
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np
import random

RADIUS = 10

PHI  = (1 + math.sqrt(5))/2

positions = [(+1, +1, +1),
            (+1, +1, -1),
            (+1, -1, +1),
            (+1, -1, -1),
            (-1, +1, +1),
            (-1, +1, -1),
            (-1, -1, +1),
            (-1, -1, -1),
            (0, +PHI, +1/PHI),
            (0, +PHI, -1/PHI),
            (0, -PHI, +1/PHI),
            (0, -PHI, -1/PHI),
            (+1/PHI, 0, +PHI),
            (+1/PHI, 0, -PHI),
            (-1/PHI, 0, +PHI),
            (-1/PHI, 0, -PHI),
            (+PHI, +1/PHI, 0),
            (+PHI, -1/PHI, 0),
            (-PHI, +1/PHI, 0),
            (-PHI, -1/PHI, 0)]

# positions = [(0, +1, +PHI),
#              (0, +1, -PHI),
#              (0, -1, +PHI),
#              (0, -1, -PHI),
#              (+1, +PHI, 0),
#              (+1, -PHI, 0),
#              (-1, +PHI, 0),
#              (-1, -PHI, 0),
#              (+PHI, 0, +1),
#              (+PHI, 0, -1),
#              (-PHI, 0, +1),
#              (-PHI, 0, -1)]



# positions = (6/math.sqrt(12)) * np.array([(+1, 0, 0),
#                                           (-1, 0, 0),
#                                           (0, +1, 0),
#                                           (0, -1, 0),
#                                           (0, 0, +1),
#                                           (0, 0, -1)])


class Vertex:
    def __init__(self, coords):
        self.x = coords[0]
        self.y = coords[1]
        self.z = coords[2]
        self.i = 0

        self.neighbors = []
        self.puzzlestates = []
        self.edges = []

        self.n = 0
        self.h = 0

    def distanceTo(self, vertex2):
        dx = self.x - vertex2.x
        dy = self.y - vertex2.y
        dz = self.z - vertex2.z

        return math.sqrt(dx**2 + dy**2 + dz**2)
    
    def vector(self):
        return np.array([self.x, self.y, self.z])
    
    def unit(self):
        return self.vector() / np.linalg.norm(self.vector())
    
    def addneighbor(self, vertex):
        self.neighbors.append(vertex)
        self.n += 1

        self.h = self.n
        if self.h % 2 != 0:
            self.h *= 2 
 
    def __repr__(self):
        return f"Vertex {self.i}"
    
class Face:
    def __init__(self, vertices):
        self.vertices = vertices
        self.edges = []
        self.body = None
        self.key = ()

        sumv = 0
        for vertex in vertices:
            sumv = sumv + vertex.vector()

        self.center = sumv / len(vertices)
        self.unitcenter = sumv / np.linalg.norm(sumv)

    def __repr__(self):
        return f"Face {self.key}"
    
class Edge:
    def __init__(self, vertex1, vertex2):
        self.vertex1 = vertex1
        self.vertex2 = vertex2
        self.faces = {'left': None, 'right': None}
        self.obstacles = {'left': False, 'right': False}
        self.key = ()

        self.midpoint = (vertex1.vector() + vertex2.vector())
        self.midpoint = self.midpoint / np.linalg.norm(self.midpoint)
        self.direction = vertex2.vector() - vertex1.vector()
        self.axis = np.cross(vertex1.vector(), vertex2.vector())

    def randomobstacles(self):
        self.obstacles = {'left': bool(random.getrandbits(1)),'right': bool(random.getrandbits(1))}

    def allobstacles(self):
        self.obstacles = {'left': True, 'right': True}

    def noobstacles(self):
        self.obstacles = {'left': False, 'right': False}

    def __repr__(self):
        return f"Edge {self.key}"

class Polyhedron:
    def __init__(self, positions):
        self.vertices = []
        self.faces = {}
        self.edges = {}
        self.graph = nx.Graph()

        i = 1
        for position in positions:
            v = Vertex(position)
            v.i = i
            self.vertices.append(v)
            self.graph.add_node(v)
            i = i + 1

        mindistance = math.inf
        for v1 in self.vertices:
            for v2 in self.vertices:
                d = v1.distanceTo(v2)
                if d < mindistance and v1 != v2:
                    mindistance = d

            for v2 in self.vertices:
                d = v1.distanceTo(v2)
                if d < mindistance + 0.00001 and v1 != v2:
                    v1.addneighbor(v2)
                    self.graph.add_edge(v1, v2)

        cycles = []
        i = 3
        while not cycles:
            cycles = list(nx.simple_cycles(self.graph, length_bound = i))
            i = i + 1

        for cycle in cycles:
            key = []
            for vertex in cycle:
                key.append(vertex.i)
            self.faces[tuple(key)] = Face(cycle)
            self.faces[tuple(key)].key = tuple(key)
 
        for edge in self.graph.edges:
            e = Edge(edge[0], edge[1])
            edge[0].edges.append(e)
            edge[1].edges.append(e)
            
            edgefaces = []
            for face in self.faces.values():
                if set(edge).issubset(set(face.vertices)):
                    edgefaces.append(face)

            v1 = np.cross(e.midpoint, edgefaces[0].unitcenter)
            v2 = np.cross(e.midpoint, edgefaces[1].unitcenter)

            if np.dot(v1, e.direction) > 0 and np.dot(v2, e.direction) < 0:
                e.faces['left'] = edgefaces[1]
                e.faces['right'] = edgefaces[0]
            elif np.dot(v1, e.direction) < 0 and np.dot(v2, e.direction) > 0:
                e.faces['left'] = edgefaces[0]
                e.faces['right'] = edgefaces[1]

            e.key = (edge[0].i, edge[1].i)
            self.edges[e.key] = e
            
    def allobstacles(self):
        for edge in self.edges.values():
            edge.allobstacles()
        
    
    def __repr__(self):
        output = ""
        for v in self.vertices:
            output = f"{output}Vertex {v.i} ({v.x}, {v.y}, {v.z}), Neighbors {v.neighbors}\n"
        return output
    
def main():
    p = Polyhedron(positions)
    print(p.graph)
    for edge in p.edges.values():
        print(edge)

    nx.draw(p.graph)
    plt.show()


    
if __name__ == "__main__":
    main()



