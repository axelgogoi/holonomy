# Written by Axel Gogoi
# https://github.com/axelgogoi/holonomy

# Creates a random holonomy puzzle within certain parameters

import math
import networkx as nx
import numpy as np
import random
import signal
import json
import os
import time
import sys


import plotly.express as px
import plotly.graph_objs as go

try:
    import polyhedron as poly
except BaseException:
    pass

# random.seed(876307)
# random.seed(324570)

# 1793192818926746468
# 893226051810471092
# 6836152261028399155
# 2627494293908858264
# 6652706253213779321
# 6550728306226543917

# seed = random.randrange(sys.maxsize)

# current map
# seed = 8450052400624775650 
# seed = 1163042937906463111

# seed = 8127425872746513940
# 718818215373588475 (long boi)
# 3019509472840850744 (25 path)
# 5226479447156815238 (2 loops)
# 1163042937906463111 (31 path)
# random.seed(seed)

# print(f"SEED: {seed}")

ALMOST0 = 0.00000001
PHI = (1 + math.sqrt(5)) / 2

positions = [(+1, +1, +1),
             (+1, +1, -1),
             (+1, -1, +1),
             (+1, -1, -1),
             (-1, +1, +1),
             (-1, +1, -1),
             (-1, -1, +1),
             (-1, -1, -1),
             (0, +PHI, +1 / PHI),
             (0, +PHI, -1 / PHI),
             (0, -PHI, +1 / PHI),
             (0, -PHI, -1 / PHI),
             (+1 / PHI, 0, +PHI),
             (+1 / PHI, 0, -PHI),
             (-1 / PHI, 0, +PHI),
             (-1 / PHI, 0, -PHI),
             (+PHI, +1 / PHI, 0),
             (+PHI, -1 / PHI, 0),
             (-PHI, +1 / PHI, 0),
             (-PHI, -1 / PHI, 0)]

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

# positions = np.array([(+1, 0, -1/math.sqrt(2)),
#                       (-1, 0, -1/math.sqrt(2)),
#                       (0, +1, +1/math.sqrt(2)),
#                       (0, -1, +1/math.sqrt(2))])

# positions = np.array([(+1, +1, +1),
#                       (+1, +1, -1),
#                       (+1, -1, +1),
#                       (+1, -1, -1),
#                       (-1, +1, +1),
#                       (-1, +1, -1),
#                       (-1, -1, +1),
#                       (-1, -1, -1)])


class PuzzleState:
    def __init__(self, vertex, heading, orientation=None):
        self.vertex = vertex
        self.heading = heading

        self.location = vertex.unit()
        self.x = self.location[0]
        self.y = self.location[1]
        self.z = self.location[2]

        if orientation is None:

            u1 = vertex.vector()
            u2 = vertex.neighbors[0].vector()
            u3 = u2 - u1

            v0 = u3 + u1 * np.dot(-u1, u3) / (np.linalg.norm(u1)**2)

            v = np.matmul(rotate3D(u1, heading * 2 * math.pi / self.vertex.h), v0)
            self.orientation = v / np.linalg.norm(v)

        else:
            self.orientation = orientation

        self.vx = self.orientation[0]
        self.vy = self.orientation[1]
        self.vz = self.orientation[2]

        self.plotlocation = self.location + 0.2 * self.orientation

    def transform(self, n):
        destinationvertex = self.vertex.neighbors[n]

        u1 = self.vertex.vector()
        u2 = destinationvertex.vector()

        axis = np.cross(u1, u2)
        axis = axis / np.linalg.norm(axis)

        theta = np.arccos(
            np.dot(u1, u2) / (np.linalg.norm(u1) * np.linalg.norm(u2)))

        v = np.matmul(rotate3D(axis, theta), self.orientation)

        # finding new heading

        u3 = destinationvertex.neighbors[0].vector()
        u4 = u3 - u2
        v0 = u4 + u2 * np.dot(-u2, u4) / (np.linalg.norm(u2)**2)


        difference = math.inf
        aligned = 0
        i = 0
        while difference > ALMOST0 and aligned <= 0:
            vi = np.matmul(rotate3D(u2, i * 2 * math.pi / destinationvertex.h), v0)
            difference = np.linalg.norm(np.cross(vi, v))
            aligned = np.dot(vi, v)
            if i < destinationvertex.h:
                i = i + 1
            else:
                break

        return PuzzleState(destinationvertex, i, v)

    def __repr__(self):
        return f"{self.vertex}-Heading {self.heading}"


class Rook:
    def __init__(self, puzzlestate, n_arms=0):
        self.puzzlestate = puzzlestate
        self.n_arms = n_arms

        self.arms = list(False for _ in range(puzzlestate.vertex.h))

        for i in range(n_arms):
            self.arms[i] = True


class HolonomyPuzzle:
    def __init__(self, polyhedron, arms = 0, seed = random.randrange(sys.maxsize)):
        self.vertices = polyhedron.vertices
        self.faces = polyhedron.faces
        self.edges = polyhedron.edges
        self.geometrygraph = polyhedron.graph
        self.puzzlegraph = nx.Graph()
        self.seed = seed

        self.rng = np.random.default_rng(seed)
            # random.seed(seed)
            # np.random.seed(seed)
            # os.environ["PYTHONHASHSEED"] = str(seed)
        print(f"SEED: {seed}")

        for vertex in self.vertices:
            for i in range(vertex.h):
                ps = PuzzleState(vertex, i)
                vertex.puzzlestates.append(ps)
                self.puzzlegraph.add_node(ps)

        for vertex in self.vertices:
            for ps in vertex.puzzlestates:
                for i in range(vertex.n):
                    psi = ps.transform(i)
                    for ps2 in self.puzzlegraph:
                        if (psi.location == ps2.location).all() and (
                                abs(psi.orientation - ps2.orientation) < ALMOST0).all():
                            if not self.puzzlegraph.has_edge(ps, ps2):
                                self.puzzlegraph.add_edge(ps, ps2)

        self.rook = Rook(self.vertices[0].puzzlestates[0], arms)
        self.obstacleslots = []

        self.start = False
        self.exit = False

        self.obstaclemap = {}
        for edge in self.edges.values():
            for side in edge.obstacles.keys():
                self.obstaclemap[(edge, side)] = set()
                self.processobstacle(edge, side)

        self.title = ""

    def addstartandexit(self, i_vertex, i_heading, f_vertex, f_neighbor):
        self.start = self.vertices[i_vertex - 1].puzzlestates[i_heading]

        self.finishvertex = self.vertices[f_vertex - 1]
        self.behindvertex = self.vertices[f_vertex - 1].neighbors[f_neighbor]

        for ps in self.finishvertex.puzzlestates:
            for ps2 in self.behindvertex.puzzlestates:
                if self.puzzlegraph.has_edge(ps, ps2):
                    self.puzzlegraph.remove_edge(ps, ps2)

        f_heading = ((f_neighbor * 2) + 2) % self.behindvertex.h
        
        self.exit = self.vertices[f_vertex - 1].puzzlestates[f_heading]



    def addobstacle(self, vertex1, vertex2, side):
        edge = self.edges[(vertex1, vertex2)]

        if edge.obstacles[side]:
            print(f"Obstacle already at {edge} on {side}")
            return False
        else:
            edge.obstacles[side] = True
            targetedges = self.obstaclemap[(edge, side)]

            skip = False
            removededges = []
            for targetedge in targetedges:
                ps1 = targetedge[0]
                ps2 = targetedge[1]
                if self.puzzlegraph.has_edge(ps1, ps2):
                    self.puzzlegraph.remove_edge(ps1, ps2)
                    removededges.append((ps1, ps2))
                    if self.start and self.exit:
                        if not nx.has_path(self.puzzlegraph, self.start, self.exit):
                            print(f"Skipping obstacle at {edge} on {side}")
                            skip = True
                            break
            
            if skip:
                for removededge in removededges:
                    self.puzzlegraph.add_edge(removededge[0], removededge[1])
            else:
                print(f"Added obstacle to {edge} on {side}")

            return True
        
    def removeobstacle(self, vertex1, vertex2, side):
        edge = self.edges[(vertex1, vertex2)]

        if not edge.obstacles[side]:
            print(f"No obstacle at {edge} on {side}")
            return False
        else:
            edge.obstacles[side] = False
            targetedges = self.obstaclemap[(edge, side)]

            for targetedge in targetedges:
                ps1 = targetedge[0]
                ps2 = targetedge[1]
                if not self.puzzlegraph.has_edge(ps1, ps2):
                    self.puzzlegraph.add_edge(ps1, ps2)

            if side == 'left':
                otherside = 'right'
            elif side == 'right':
                otherside = 'left'

            if edge.obstacles[otherside]:
                replacementedges = self.obstaclemap[(edge, otherside)]
                for replacementedge in replacementedges:
                    ps1 = replacementedge[0]
                    ps2 = replacementedge[1]
                    if self.puzzlegraph.has_edge(ps1, ps2):
                        self.puzzlegraph.remove_edge(ps1, ps2)

            for ps in self.finishvertex.puzzlestates:
                for ps2 in self.behindvertex.puzzlestates:
                    if self.puzzlegraph.has_edge(ps, ps2):
                        self.puzzlegraph.remove_edge(ps, ps2)

            print(f"Removed obstacle from {edge} on {side}")

            return True

    def addrandomobstacle(self):
        edge = self.rng.choice(list(self.edges.values()))
        side = self.rng.choice(['left', 'right'])

        result = self.addobstacle(edge.key[0], edge.key[1], side)
        return (result, edge, side)
    
    def removerandomobstacle(self):
        edge = self.rng.choice(list(self.edges.values()))
        side = self.rng.choice(['left', 'right'])

        result = self.removeobstacle(edge.key[0], edge.key[1], side)
        return (result, edge, side)
    
    def processobstacle(self, edge, side):
        result = True

        for ps in edge.vertex1.puzzlestates:
            x1 = edge.vertex1.unit()
            x2 = edge.vertex2.unit()
            v1 = ps.orientation

            u1 = np.cross(np.cross(x1, v1), np.cross(x2, v1))
            u2 = np.dot(np.cross(x1, x2), np.cross(x1, v1))

            if np.linalg.norm(u1) < ALMOST0 and u2 > 0:
                h0 = ps.heading

        ps1edge = {}
        for ps in edge.vertex1.puzzlestates:
            for ps2 in edge.vertex2.puzzlestates:
                if self.puzzlegraph.has_edge(ps, ps2):
                    ps1edge[ps.heading] = (ps, ps2)

        h = edge.vertex1.h
        for i in range(len(self.rook.arms)):
            if self.rook.arms[i]:
                if side == 'left':
                    leftindices = np.array(range(h0 + 1, int(h0 + h/2))) % h
                    for l in leftindices:
                        if (l - i)%h in ps1edge.keys():
                            ps1 = ps1edge[(l - i)%h][0]
                            ps2 = ps1edge[(l - i)%h][1]
                            self.obstaclemap[(edge, side)].add((ps1, ps2))

                if side == 'right':
                    rightindices = np.array(range(int(h0 + 1 - h/2), h0)) % h
                    for r in rightindices:
                        if (r - i)%h in ps1edge.keys():
                            ps1 = ps1edge[(r - i)%h][0]
                            ps2 = ps1edge[(r - i)%h][1]
                            self.obstaclemap[(edge, side)].add((ps1, ps2))

        return result
    
    def populate(self, pruning_enabled = True, prune_attempts = 100, pathrange = (1, 1)):
        if self.start and self.exit:
            print("populating...")

            self.clearobstacles()
            n_paths = math.inf
            while n_paths > pathrange[1]:
                self.addrandomobstacle()
                with timeout(seconds = 0.02, error_message="Timed out, too many paths"):
                    try:
                        paths = nx.all_simple_paths(self.puzzlegraph, self.start, self.exit)
                        n_paths = len(list(paths))
                        print(f"Number of paths: {n_paths}")
                    except Exception as ex:
                        print(ex)
                        pass

            if pruning_enabled:
                print("pruning...")
                report = ""
                count = 0
                for _ in range(prune_attempts):
                    (result, edge, side) = self.removerandomobstacle()
                    failed = False
                    if result:
                        with timeout(seconds = 0.2, error_message="Timed out, restarting population"):
                            try:
                                paths = nx.all_simple_paths(self.puzzlegraph, self.start, self.exit)
                                n_paths = len(list(paths))
                                if not pathrange[0] <= n_paths <= pathrange[1]:
                                    self.addobstacle(edge.key[0], edge.key[1], side)
                                else:
                                    count += 1
                                    report += f"Pruned obstacle from {edge}, on {side}\n"
                            except Exception as ex:
                                print(ex)
                                failed = True
                    if failed:
                        self.addobstacle(edge.key[0], edge.key[1], side)
                        break

                print("\n")
                print(report)
                return report
            else:
                return "Nothing to report"

        else:
            print("No start or exit")

    def addrandomobstacles(self, attempts):
        if self.start and self.exit:
            for _ in range(attempts):
                self.addrandomobstacle()
        else:
            print("No start or exit")

    def getobstaclecount(self):
        output = ''
        count = 0
        for edge in self.edges.values():
            for side in edge.obstacles.keys():
                if edge.obstacles[side]:
                    count += 1
                    output += f"Obstacle at {edge}, on {side}\n"

        output += f"Total number of obstacles: {count}\n"
        print(output)
        return output
    
    def getpath(self, display = True):
        if self.start and self.exit:
            with timeout(seconds = 0.2, error_message="Timed out, too many paths"):
                try:
                    paths = list(nx.all_simple_paths(self.puzzlegraph, self.start, self.exit))
                    if paths:
                        paths.sort(key = len)
                        if display:
                            for i in range(len(paths)):
                                print(f"Path {i}, length {len(paths[i])}")
                                for step in paths[i]:
                                    print(step)
                                print(f"Path branchiness: {self.pathbranchiness(paths[i])}")
                        print(f"{len(paths)} paths between start and exit")
                    else:
                        print("No path exists between start and exit")
                    return paths
                
                except Exception as ex:
                    print(ex)
                    return True

    def getcustompath(self, i_vertex, i_heading, f_vertex, f_heading):
        start = self.vertices[i_vertex-1].puzzlestates[i_heading]
        finish = self.vertices[f_vertex-1].puzzlestates[f_heading]

        haspath = nx.has_path(self.puzzlegraph, start, finish)
        if haspath:
            paths = nx.all_simple_paths(self.puzzlegraph, start, finish)
            paths = list(paths)
            for i in range(len(paths)):
                print(f"Path {i}, length {len(paths[i])}")
                for step in paths[i]:
                    print(step)
            path = paths[0]
            return path
        else:
            print("No path exists between start and exit")
    
    def moverook(self, move):
        currentvertex = self.rook.puzzlestate.vertex
        targetvertex = self.rook.puzzlestate.vertex.neighbors[move]
        puzzleedges = self.puzzlegraph.edges(self.rook.puzzlestate)

        targetstate = None
        for puzzleedge in puzzleedges:
            for ps in puzzleedge:
                if ps.vertex == targetvertex:
                    targetstate = ps

        if targetstate is not None:
            self.rook.puzzlestate = targetstate
            print(f"Rook moved from {currentvertex} to {targetvertex}")
        else:
            print("move not allowed")

    def saveobstacles(self, filename = 'obstacles.json'):
        saved_obstacles = []
        for edgekey in self.edges.keys():
            for side in self.edges[edgekey].obstacles.keys():
                if self.edges[edgekey].obstacles[side]:
                    obstacle = {}
                    obstacle['v1'] = edgekey[0]
                    obstacle['v2'] = edgekey[1]
                    obstacle['side'] = side
                    saved_obstacles.append(obstacle)
        
        self.title = filename
        with open(filename, 'w') as f:
            json.dump(saved_obstacles, f)  

    def loadobstacles(self, filename = 'obstacles.json'):

        with open(filename, 'r') as f:
            self.clearobstacles()
            saved_obstacles = json.load(f)
            
            print("\nLoading obstacles...")
            for o in saved_obstacles:
                self.addobstacle(o['v1'], o['v2'], o['side'])
            print("")
            self.title = filename
            return self.getobstaclecount()
        
    def clearobstacles(self):
        print("Clearing obstacles...")
        for edgekey in self.edges.keys():
            for side in self.edges[edgekey].obstacles.keys():
                if self.edges[edgekey].obstacles[side]:
                    self.removeobstacle(edgekey[0], edgekey[1], side)

    def maingraph(self):
        components = sorted(nx.connected_components(self.puzzlegraph), key=len, reverse=True)
        maingraph = self.puzzlegraph.subgraph(components[0])
        
        return maingraph
    
    def cycles(self):
        with timeout(seconds = 0.1, error_message="Timed out, too many loops"):
            try:
                n_cycles = len(list(nx.simple_cycles(self.puzzlegraph)))
                return n_cycles
            except Exception as ex:
                print(ex)
                return -1
            
            
    def exitfactor(self):
        maingraph = self.maingraph()
        d = nx.diameter(maingraph)
        r = nx.radius(maingraph)
        e = nx.eccentricity(maingraph, self.exit)

        return (e - r) / (d - r)
    
    def branchiness(self):
        maingraph = self.maingraph()
        branchiness = 0
        for node in self.maingraph():
            if len(list(maingraph.neighbors(node))) == 1:
                branchiness += 1
            elif len(list(maingraph.neighbors(node))) == 2:
                branchiness -= 1
            else:
                branchiness += len(list(maingraph.neighbors(node))) - 1
        return branchiness
    
    def pathbranchiness(self, path):
        n = len(path)
        maingraph = self.maingraph()
        branches = 0
        for node in path:
            if len(list(maingraph.neighbors(node))) > 2:
                branches += 1

        return branches / n
            

    def info(self):
        n_cycles = self.cycles()
        maingraph = self.maingraph()

        print(f"\nSEED: {self.seed}")
        print(f"Complete graph: {self.puzzlegraph}")
        print(f"Connected graph: {maingraph}")
        print(f"Number of cycles: {n_cycles}")
        print(f"Diameter: {nx.diameter(maingraph)}")
        print(f"Radius: {nx.radius(maingraph)}")
        print(f"Exit eccentricity: {nx.eccentricity(maingraph, self.exit)}")
        print(f"Exit factor: {self.exitfactor()}")
        print(f"Branchiness: {self.branchiness()}\n")

    def search(self, pathrange, cyclerange, minlength, exitrange, minbranch, p_e = True, minsize = 105, maxdiam = 35): 
        attempts = 0
        while True:
            attempts += 1
            self.populate(pruning_enabled = p_e, pathrange = pathrange)

            paths = self.getpath()
            pause = 0.15
            self.getobstaclecount()
            maingraph = self.maingraph()

        
            if not paths or type(paths) is not list:
                print(f"Attempt {attempts} failed, no path or too many paths")
                time.sleep(pause)

            elif len(maingraph) < minsize:
                print(f"Attempt {attempts} failed, graph is too small")
                time.sleep(pause)

            elif nx.diameter(maingraph) > maxdiam:
                print(f"Attempt {attempts} failed, diameter is too large")
                time.sleep(pause)

            elif len(list(maingraph.neighbors(self.exit))) > 1:
                print(f"Attempt {attempts} failed, wasted nodes around exit")
                time.sleep(pause)

            elif not pathrange[0] <= len(paths) <= pathrange[1]:
                print(f"Attempt {attempts} failed, paths not in range {pathrange}")
                time.sleep(pause)

            elif not cyclerange[0] <= self.cycles() <= cyclerange[1]:
                print(f"Attempt {attempts} failed, {self.cycles()} loops not in range {cyclerange}")
                time.sleep(pause)

            elif len(paths[0]) < minlength:
                print(f"Attempt {attempts} failed, start and exit are too close")
                time.sleep(pause)
                
            elif not exitrange[0] < self.exitfactor() < exitrange[1]:
                print(f"Attempt {attempts} failed, exit is not central enough")
                time.sleep(pause)
            
            elif self.branchiness() < minbranch:
                print(f"Attempt {attempts} failed, graph does not have enough branches")
                time.sleep(pause)

            elif self.pathbranchiness(paths[0]) < 0.45:
                print(f"Attempt {attempts} failed, path does not have enough branches")
                time.sleep(pause)

            else:
                break

        print(f"Solution found after {attempts} attempts")

class timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.setitimer(signal.ITIMER_REAL, self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)


def rotate3D(axis, theta):
    axis = axis / (np.linalg.norm(axis))
    a = math.cos(theta / 2.0)
    b, c, d = - axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def main():
    p = poly.Polyhedron(positions)
    h = HolonomyPuzzle(p, 2)
    h.addstartandexit(1, 0, 1, 0)

    # print(h.puzzlegraph)

    # h.search((1, 3), (1, 7), 22, (0.5, 0.75), 37)
    h.loadobstacles(filename='obstacles37.json')
 
    h.getobstaclecount()
    paths = h.getpath()
    h.info()
    
    # h.saveobstacles(filename = 'obstacles22.json')

    plotmap(h, paths[0])
    plotpuzzle(h, paths[0])
    

def plotmap(h, path = None):
    points = nx.kamada_kawai_layout(h.maingraph(), dim = 3)
    annotations = []
    colors = []
    x, y, z = [], [], []
    for ps in points.keys():
        x.append(points[ps][0])
        y.append(points[ps][1])
        z.append(points[ps][2])

        u = ps.vertex.unit()
        i = (np.arctan2(u[0],u[2]) + np.pi)/(2*np.pi)
        color = px.colors.sample_colorscale(px.colors.cyclical.HSV, i)

        annotation = dict(font_color = color[0], font_size = 11, font_family="Arial Black", showarrow=False)
        annotation['x'] = points[ps][0]
        annotation['y'] = points[ps][1]
        annotation['z'] = points[ps][2]
        annotation['text'] = f"{ps.vertex.i},{ps.heading}"
        annotation['xanchor'] = 'left'
        annotation['yanchor'] = 'bottom'
        annotations.append(annotation)

        colors.append(color[0])

    fig_map = px.scatter_3d(x=x, y=y, z=z, template="plotly_dark", title = h.title)
    fig_map.update_traces(marker_size = 5)
    fig_map.update_traces(marker_color = colors)

    x = [points[h.start][0]]
    y = [points[h.start][1]]
    z = [points[h.start][2]]
    fig_map.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color = 'lime', size = 12, opacity = 0.6)))

    x = [points[h.exit][0]]
    y = [points[h.exit][1]]
    z = [points[h.exit][2]]
    fig_map.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color = 'red', size = 12, opacity = 0.6)))
    
    for edge in h.puzzlegraph.edges:
        if edge[0] in h.maingraph() and edge[1] in h.maingraph():
            x, y, z = [], [], []

            x = [points[edge[0]][0], points[edge[1]][0]]
            y = [points[edge[0]][1], points[edge[1]][1]]
            z = [points[edge[0]][2], points[edge[1]][2]]
            
            line=dict(width=4, color='#888888')
            fig_map.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', showlegend=False, line=line))

    if path:
        x, y, z = [], [], []
        for ps in path:
            x.append(points[ps][0])
            y.append(points[ps][1])
            z.append(points[ps][2])

        line = dict(width = 5, color = 'white')
        fig_map.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', showlegend=False, line=line))

    annotations.reverse()
    fig_map.update_layout(scene=dict(annotations=annotations))
    fig_map.show()


def plotpuzzle(h, path = None):
    x, y, z = [], [], []
    annotations = []
    for ps in h.puzzlegraph:
        x.append(ps.plotlocation[0])
        y.append(ps.plotlocation[1])
        z.append(ps.plotlocation[2])

        annotation = dict(font_color = "black", font_size = 11, font_family="Arial Black", showarrow=False)
        annotation['x'] = ps.plotlocation[2]
        annotation['y'] = ps.plotlocation[0]
        annotation['z'] = ps.plotlocation[1]
        annotation['text'] = ps.heading
        annotations.append(annotation)

    fig_puzzle = px.scatter_3d(x=z,y=x,z=y,color=np.arctan2(x, z),color_continuous_scale=px.colors.cyclical.HSV,template="plotly_dark", title = h.title)
    fig_puzzle.update_traces(marker_size = 10)

    x, y, z = [], [], []

    if path:
        for ps in path:
            x.append(ps.plotlocation[0])
            y.append(ps.plotlocation[1])
            z.append(ps.plotlocation[2])

        line = dict(width=20, color=np.arctan2(x, z), cmin=-np.pi, cmax=np.pi, colorscale=px.colors.cyclical.HSV)
        fig_puzzle.add_trace(go.Scatter3d(x=z,y=x,z=y,mode='lines',showlegend=False, line = line))

    x, y, z = [], [], []
    for edge in h.edges.values():
        if edge.obstacles['left'] or edge.obstacles['right']:
            midpoint = edge.vertex1.vector() + edge.vertex2.vector()
            midpoint = midpoint / np.linalg.norm(midpoint)
            direction = np.cross(
                edge.vertex1.vector(),
                edge.vertex2.vector())

            if edge.obstacles['left']:
                v = midpoint + 0.06 * direction
                v = v / np.linalg.norm(v)
                x.append(v[0])
                y.append(v[1])
                z.append(v[2])
            if edge.obstacles['right']:
                v = midpoint - 0.06 * direction
                v = v / np.linalg.norm(v)
                x.append(v[0])
                y.append(v[1])
                z.append(v[2])

    marker = dict(color='white', size=15)
    fig_puzzle.add_trace(go.Scatter3d(x=z,y=x,z=y,mode='markers',showlegend=False, marker = marker))

    for edge in h.puzzlegraph.edges:
        x, y, z = [], [], []

        x = [edge[0].plotlocation[0], edge[1].plotlocation[0]]
        y = [edge[0].plotlocation[1], edge[1].plotlocation[1]]
        z = [edge[0].plotlocation[2], edge[1].plotlocation[2]]

        line = dict(width=5, color=np.arctan2(x, z), cmin=-np.pi, cmax=np.pi, colorscale=px.colors.cyclical.HSV)
        fig_puzzle.add_trace(go.Scatter3d(x=z,y=x,z=y,mode='lines',showlegend=False,line = line))     
   
    for vertex in h.vertices:
        u = vertex.unit() * 1.2
        annotations.append(dict(font_size = 20, font_family="Avenir Next", showarrow=False, x = u[2], y = u[0], z = u[1], text = vertex.i))

    annotations.reverse()
    fig_puzzle.update_layout(scene=dict(annotations=annotations))
    fig_puzzle.show()


if __name__ == "__main__":
    main()
