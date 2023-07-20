# Written by Axel Gogoi
# https://github.com/axelgogoi/holonomy

# Models a holonomy puzzle in Fusion 360 from a stored list of obstacles

import adsk.core, adsk.fusion, adsk.cam, traceback
import math
import time
import os, sys
from statistics import median

path = os.path.dirname(os.path.realpath("/usr/local/lib/python3.9/site-packages/networkx"))

if not path in sys.path:
   sys.path.append(path) 

import numpy as np
from . import polyhedron as poly
from . import holonomy as holo

RADIUS = 7
RATIO = 0.68
FLOORTHICKNESS = 0.2
OUTERHOLE = 0.82
INNERHOLE = 0.85

CUSHIONGAP = 0.13 * 2.54
CORETHICKNESS = 0.26
SOCKETDEPTH = 0.13
SOCKETSIZE = 1.79 + 0.017
SOCKETTHICKNESS = 0.3

GROOVEFACTOR = 1.05

INNERRADIUS = RADIUS * RATIO

GROOVEDEPTH = 0.8
LIPHEIGHT = 0.13
RIDGEHEIGHT = 0.55
OBSTACLEHEIGHT = 0.35

ZERO = adsk.core.Point3D.create(0, 0, 0)    # defining zero point for convenience

PHI  = (1 + math.sqrt(5))/2

positions = ((math.sqrt(3 - PHI))) / PHI * np.array([(+1, +1, +1),
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
                                                     (-PHI, -1/PHI, 0)])

HOLLOWCORE = True
GROOVES    = True
OBSTACLES  = True
POCKETS    = True
  

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui  = app.userInterface

        design = adsk.fusion.Design.cast(app.activeProduct)
        occs = design.rootComponent.occurrences
        mat = adsk.core.Matrix3D.create()
        newOcc = occs.addNewComponent(mat)        
        comp = adsk.fusion.Component.cast(newOcc.component)

        sketches = comp.sketches
        revolves = comp.features.revolveFeatures
        extrudes = comp.features.extrudeFeatures
        lofts    = comp.features.loftFeatures
        offsets  = comp.features.offsetFeatures
        thickens = comp.features.thickenFeatures
        sweeps   = comp.features.sweepFeatures
        mirrors  = comp.features.mirrorFeatures
        cplanes  = comp.constructionPlanes
        cpoints  = comp.constructionPoints
        combines = comp.features.combineFeatures

        timelineGroups = design.timeline.timelineGroups
        
        xyPlane = comp.xYConstructionPlane
        sphere_sk = sketches.add(xyPlane)     # creating sphere sketch
        sphere_sk.sketchCurves.sketchArcs.addByThreePoints(point2d(0, RADIUS), point2d(RADIUS, 0), point2d(0, -RADIUS))
        axis = sphere_sk.sketchCurves.sketchLines.addByTwoPoints(point2d(0,RADIUS), point2d(0,-RADIUS))
        centerpoint = sphere_sk.sketchPoints.add(ZERO)

        sphereprof = adsk.core.ObjectCollection.create()
        for profile in sphere_sk.profiles:
            sphereprof.add(profile)


        dodecahedron = poly.Polyhedron(positions)
        h = holo.HolonomyPuzzle(dodecahedron, 2)
        h.addstartandexit(1, 0, 1, 0)
        report = h.loadobstacles("/Users/axelgogoi/Library/Application Support/Autodesk/Autodesk Fusion 360/API/Scripts/holonomy/obstacles39.json")
        app.log(report)


        facebodies = []
        corebodies = []
        hollowbodies = []
        rays_sk = sketches.add(xyPlane)

        faces = h.faces.values()

        progress1 = Progress(ui, "Progress - Building Elements", faces)

        for face in faces:
            progress1.check()
            
            revolveIn = revolves.createInput(sphereprof, axis, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
            angle = adsk.core.ValueInput.createByReal(2*math.pi)
            revolveIn.setAngleExtent(False, angle)
            sphere = revolves.add(revolveIn)

            face_sk = sketches.add(xyPlane)
            n = len(face.vertices)
            for i in range(n):
                v1 = 1.1 * RADIUS * face.vertices[i].vector()
                v2 = 1.1 * RADIUS * face.vertices[(i+1)%n].vector()
                face_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(v1), point3d(v2))   

                v1 = INNERRADIUS * face.vertices[i].vector()
                v2 = INNERRADIUS * face.vertices[(i+1)%n].vector()
                face_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(v1), point3d(v2))

                v1 = (INNERRADIUS - CUSHIONGAP) * face.vertices[i].vector()
                v2 = (INNERRADIUS - CUSHIONGAP) * face.vertices[(i+1)%n].vector()
                face_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(v1), point3d(v2))

                v1 = (INNERRADIUS - CUSHIONGAP - CORETHICKNESS - SOCKETDEPTH) * face.vertices[i].vector()
                v2 = (INNERRADIUS - CUSHIONGAP - CORETHICKNESS - SOCKETDEPTH) * face.vertices[(i+1)%n].vector()
                face_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(v1), point3d(v2))

                u1 = (face.vertices[i].vector() - face.center) 
                u2 = (face.vertices[(i+1)%n].vector() - face.center)

                v0 = (INNERRADIUS - CUSHIONGAP - CORETHICKNESS - SOCKETDEPTH) * face.center

                v1 = v0 + SOCKETSIZE * u1 / np.linalg.norm(u1)
                v2 = v0 + SOCKETSIZE * u2 / np.linalg.norm(u2)
                face_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(v1), point3d(v2))

                # v1 = v0 + (SOCKETSIZE + SOCKETTHICKNESS) * u1 / np.linalg.norm(u1)
                # v2 = v0 + (SOCKETSIZE + SOCKETTHICKNESS) * u2 / np.linalg.norm(u2)
                # face_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(v1), point3d(v2))

            face_profiles = sorted(face_sk.profiles, key=lambda profile: profile.areaProperties().area, reverse=True)
            
            outerprof = adsk.core.ObjectCollection.create()
            innerprof = adsk.core.ObjectCollection.create()
            coreprof  = adsk.core.ObjectCollection.create()
            
            outerprof.add(face_profiles[0])
            innerprof.add(face_profiles[1])
            coreprof.add(face_profiles[2])

            if not corebodies:
                loftIn = lofts.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
            else:
                loftIn = lofts.createInput(adsk.fusion.FeatureOperations.JoinFeatureOperation)
                loftIn.participantBodies = corebodies

            loftIn.loftSections.add(coreprof)
            loftIn.loftSections.add(centerpoint)
            loft = lofts.add(loftIn)
            corebodies.append(loft.bodies.item(0))

            if HOLLOWCORE:
                if not hollowbodies:
                    loftIn = lofts.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
                else:
                    loftIn = lofts.createInput(adsk.fusion.FeatureOperations.JoinFeatureOperation)
                    loftIn.participantBodies = hollowbodies

                hollowprof = adsk.core.ObjectCollection.create()
                socketprof = adsk.core.ObjectCollection.create()
                hollowprof.add(face_profiles[3])
                socketprof.add(face_profiles[-1])

                loftIn.loftSections.add(hollowprof)
                loftIn.loftSections.add(centerpoint)
                loft = lofts.add(loftIn)
                hollowbodies.append(loft.bodies.item(0))

                extrudeIn = extrudes.createInput(socketprof, adsk.fusion.FeatureOperations.JoinFeatureOperation)
                distance = adsk.core.ValueInput.createByReal(SOCKETDEPTH)
                extrudeIn.setSymmetricExtent(distance, False)
                extrudeIn.participantBodies = hollowbodies
                extrudes.add(extrudeIn)

            loftIn = lofts.createInput(adsk.fusion.FeatureOperations.IntersectFeatureOperation)
            loftIn.loftSections.add(innerprof)
            loftIn.loftSections.add(outerprof)
            loftIn.participantBodies = [sphere.bodies.item(0)]
            loft = lofts.add(loftIn)
            face.body = loft.bodies.item(0)
            facebodies.append(face.body)

            progress1.update()

        progress1.finish()

        if OBSTACLES or GROOVES:

            edges_sk = sketches.add(xyPlane)
            edges = h.edges.values()
            progress2 = Progress(ui,"Progress - Processing Edges", edges)

            for edge in edges:
                progress2.check()
                
                ### OBSTACLES

                if OBSTACLES:

                    if edge.obstacles['left'] or edge.obstacles['right']:
                        obstacle_sk = sketches.add(xyPlane)

                        a = RADIUS * edge.midpoint
                        b1 = np.matmul(holo.rotate3D(edge.axis, +5 * np.pi/180), a)
                        b2 = np.matmul(holo.rotate3D(edge.axis, -5 * np.pi/180), a)

                        c  = (1 + OBSTACLEHEIGHT/RADIUS) * a
                        d1 = np.matmul(holo.rotate3D(edge.axis, +4 * np.pi/180), c)
                        d2 = np.matmul(holo.rotate3D(edge.axis, -4 * np.pi/180), c)

                        obstacle_sk.sketchCurves.sketchArcs.addByThreePoints(point3d(b1), point3d(a), point3d(b2))
                        obstacle_sk.sketchCurves.sketchArcs.addByThreePoints(point3d(d1), point3d(c), point3d(d2))
                        obstacle_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(b1), point3d(d1))
                        obstacle_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(b2), point3d(d2))

                        l1 = RADIUS * edge.faces['left'].unitcenter
                        l2 = RADIUS * (l1 + a) / np.linalg.norm(l1 + a)

                        r1 = RADIUS * edge.faces['right'].unitcenter
                        r2 = RADIUS * (l1 + a) / np.linalg.norm(l1 + a)

                        if edge.obstacles['left']:
                            l1 = RADIUS * edge.faces['left'].unitcenter
                            l2 = RADIUS * (l1 + a) / np.linalg.norm(l1 + a)

                            leftpath = obstacle_sk.sketchCurves.sketchArcs.addByThreePoints(point3d(l1), point3d(l2), point3d(a))
                            leftpath = comp.features.createPath(leftpath, False)

                            obstacleprof = adsk.core.ObjectCollection.create()
                            for profile in obstacle_sk.profiles:
                                obstacleprof.add(profile)

                            sweepIn = sweeps.createInput(obstacleprof, leftpath, adsk.fusion.FeatureOperations.JoinFeatureOperation)
                            sweepIn.participantBodies = [edge.faces['left'].body]
                            sweeps.add(sweepIn)

                        if edge.obstacles['right']:
                            r1 = RADIUS * edge.faces['right'].unitcenter
                            r2 = RADIUS * (r1 + a) / np.linalg.norm(r1 + a)

                            rightpath = obstacle_sk.sketchCurves.sketchArcs.addByThreePoints(point3d(r1), point3d(r2), point3d(a))
                            rightpath = comp.features.createPath(rightpath, False)

                            obstacleprof = adsk.core.ObjectCollection.create()
                            for profile in obstacle_sk.profiles:
                                obstacleprof.add(profile)

                            sweepIn = sweeps.createInput(obstacleprof, rightpath, adsk.fusion.FeatureOperations.JoinFeatureOperation)
                            sweepIn.participantBodies = [edge.faces['right'].body]
                            sweeps.add(sweepIn)

                ### GROOVES

                if GROOVES:

                    arc1 = point3d(RADIUS * edge.vertex1.unit())
                    arc2 = point3d(RADIUS * edge.midpoint)
                    arc3 = point3d(RADIUS * edge.vertex2.unit())

                    arc = edges_sk.sketchCurves.sketchArcs.addByThreePoints(arc1, arc2, arc3)

                    cplaneIn = cplanes.createInput()
                    position = adsk.core.ValueInput.createByReal(0)
                    cplaneIn.setByDistanceOnPath(arc, position)
                    cplane = cplanes.add(cplaneIn)

                    groove_sk = sketches.add(cplane)
                    projected = groove_sk.project(centerpoint)
                    projalignpoint = projected.item(0).geometry.asArray()

                    a1 = point2d(r = RADIUS - GROOVEDEPTH, theta = +1.5 * GROOVEFACTOR, origin = projalignpoint)
                    b  = point2d(r = RADIUS - GROOVEDEPTH, theta = 0, origin = projalignpoint)
                    a2 = point2d(r = RADIUS - GROOVEDEPTH, theta = -1.5 * GROOVEFACTOR, origin = projalignpoint)

                    c1 = point2d(r = RADIUS - GROOVEDEPTH +  RIDGEHEIGHT, theta = +2 * GROOVEFACTOR, origin = projalignpoint)
                    c2 = point2d(r = RADIUS - GROOVEDEPTH +  RIDGEHEIGHT, theta = -2 * GROOVEFACTOR, origin = projalignpoint)

                    d1 = point2d(r = RADIUS - GROOVEDEPTH +  RIDGEHEIGHT, theta = +2.5 * GROOVEFACTOR, origin = projalignpoint)
                    d2 = point2d(r = RADIUS - GROOVEDEPTH +  RIDGEHEIGHT, theta = -2.5 * GROOVEFACTOR, origin = projalignpoint)

                    e1 = point2d(r = RADIUS - GROOVEDEPTH +  RIDGEHEIGHT, theta = +3.5 * GROOVEFACTOR, origin = projalignpoint)
                    e2 = point2d(r = RADIUS - GROOVEDEPTH +  RIDGEHEIGHT, theta = -3.5 * GROOVEFACTOR, origin = projalignpoint)

                    f1 = point2d(r = RADIUS - GROOVEDEPTH, theta = +7 * GROOVEFACTOR, origin = projalignpoint)
                    f2 = point2d(r = RADIUS - GROOVEDEPTH, theta = -7 * GROOVEFACTOR, origin = projalignpoint)

                    g1 = point2d(r = RADIUS - GROOVEDEPTH, theta = +8 *  GROOVEFACTOR, origin = projalignpoint)
                    g2 = point2d(r = RADIUS - GROOVEDEPTH, theta = -8 * GROOVEFACTOR, origin = projalignpoint)

                    h1 = point2d(r = RADIUS - GROOVEDEPTH, theta = +11 * GROOVEFACTOR, origin = projalignpoint)
                    h2 = point2d(r = RADIUS - GROOVEDEPTH, theta = -11 * GROOVEFACTOR, origin = projalignpoint)

                    i1 = point2d(r = RADIUS - LIPHEIGHT, theta = +7 * GROOVEFACTOR, origin = projalignpoint)
                    i2 = point2d(r = RADIUS - LIPHEIGHT, theta = -7 * GROOVEFACTOR, origin = projalignpoint)

                    k1 = point2d(r = 1.25 * RADIUS, theta = +13 * GROOVEFACTOR, origin = projalignpoint)
                    k2 = point2d(r = 1.25 * RADIUS, theta = -13 * GROOVEFACTOR, origin = projalignpoint)

                    r1 = point2d(r = 0, theta = 0, origin = projalignpoint)
                    r2 = point2d(x = 0, y = RADIUS, origin = projalignpoint)

                    groove_sk.sketchCurves.sketchArcs.addByThreePoints(a1, b, a2)

                    groove_sk.sketchCurves.sketchLines.addByTwoPoints(a1, c1)
                    groove_sk.sketchCurves.sketchArcs.addByThreePoints(c1, d1, e1)
                    groove_sk.sketchCurves.sketchLines.addByTwoPoints(e1, f1)
                    groove_sk.sketchCurves.sketchArcs.addByThreePoints(f1, g1, h1)
                    groove_sk.sketchCurves.sketchLines.addByTwoPoints(h1, i1)
                    groove_sk.sketchCurves.sketchLines.addByTwoPoints(i1, k1)

                    groove_sk.sketchCurves.sketchLines.addByTwoPoints(a2, c2)
                    groove_sk.sketchCurves.sketchArcs.addByThreePoints(c2, d2, e2)
                    groove_sk.sketchCurves.sketchLines.addByTwoPoints(e2, f2)
                    groove_sk.sketchCurves.sketchArcs.addByThreePoints(f2, g2, h2)
                    groove_sk.sketchCurves.sketchLines.addByTwoPoints(h2, i2)
                    groove_sk.sketchCurves.sketchLines.addByTwoPoints(i2, k2)

                    groove_sk.sketchCurves.sketchLines.addByTwoPoints(k1, k2)

                    axis = groove_sk.sketchCurves.sketchLines.addByTwoPoints(r1, r2)

                    grooveprof = adsk.core.ObjectCollection.create()
                    for profile in groove_sk.profiles:
                        grooveprof.add(profile)

                    revolveIn = revolves.createInput(grooveprof, axis, adsk.fusion.FeatureOperations.CutFeatureOperation)
                    angle = adsk.core.ValueInput.createByReal(2*math.pi)
                    revolveIn.setAngleExtent(False, angle)
                    revolveIn.participantBodies = [edge.faces['left'].body, edge.faces['right'].body]
                    revolves.add(revolveIn)

                progress2.update()

            progress2.finish()

            vertices = h.vertices
            progress3 = Progress(ui,"Progress - Cleaning Up Vertices", vertices)
            for vertex in vertices:
                progress3.check()

                cleanup_sk = sketches.add(xyPlane)

                a = vertex.unit()
                axis = vertex.edges[0].axis
                c = np.matmul(holo.rotate3D(axis, +9.3 * np.pi/180), a)
                e = np.matmul(holo.rotate3D(axis, +7.3 * np.pi/180), a)

                b = (a + c) / np.linalg.norm(a + c)
                d = (a + e) / np.linalg.norm(a + e)

                a1 = point3d((RADIUS - GROOVEDEPTH) * a)
                b = point3d((RADIUS - GROOVEDEPTH) * b)
                c = point3d((RADIUS - GROOVEDEPTH) * c)

                a2 = point3d(RADIUS * a)
                d = point3d(RADIUS * d)
                e = point3d(RADIUS * e)

                axis2 = cleanup_sk.sketchCurves.sketchLines.addByTwoPoints(a1, a2)
                cleanup_sk.sketchCurves.sketchLines.addByTwoPoints(c, e)
                cleanup_sk.sketchCurves.sketchArcs.addByThreePoints(a1, b, c)
                cleanup_sk.sketchCurves.sketchArcs.addByThreePoints(a2, d, e)

                cleanupprof = adsk.core.ObjectCollection.create()
                for profile in cleanup_sk.profiles:
                    cleanupprof.add(profile)

                revolveIn = revolves.createInput(cleanupprof, axis2, adsk.fusion.FeatureOperations.CutFeatureOperation)
                angle = adsk.core.ValueInput.createByReal(2*math.pi)
                revolveIn.setAngleExtent(False, angle)
                sphere = revolves.add(revolveIn)

                progress3.update()

            progress3.finish()


        ## POCKETS AND HOLES

        if POCKETS:
            progress4 = Progress(ui,"Progress - Processing Faces", faces) 
            for face in faces:
                progress4.check()

                pocket_sk = sketches.add(xyPlane)
                k = 0.44

                n = len(face.vertices)
                for i in range(n):
                    v1 = (INNERRADIUS + FLOORTHICKNESS) * ((1-k)*face.unitcenter + k * face.vertices[i].vector())
                    v2 = (INNERRADIUS + FLOORTHICKNESS) * ((1-k)*face.unitcenter + k * face.vertices[(i+1)%n].vector())
                   
                    pocket_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(v1),  point3d(v2))

                for i in range(n):
                    v1 = (INNERRADIUS + FLOORTHICKNESS) * ((1-k)*face.unitcenter + k * 4 * face.vertices[i].vector())
                    v2 = (INNERRADIUS + FLOORTHICKNESS) * ((1-k)*face.unitcenter + k * 4 * face.vertices[(i+1)%n].vector())

                    pocket_sk.sketchCurves.sketchLines.addByTwoPoints(point3d(v1),  point3d(v2))

                innerprof = adsk.core.ObjectCollection.create()
                outerprof = adsk.core.ObjectCollection.create()

                innerprof.add(pocket_sk.profiles[0])
                outerprof.add(pocket_sk.profiles[1])

                loftIn = lofts.createInput(adsk.fusion.FeatureOperations.CutFeatureOperation)
                loftIn.loftSections.add(innerprof)
                loftIn.loftSections.add(outerprof)
                loftIn.participantBodies = [face.body]
                loft = lofts.add(loftIn)
                face.body = loft.bodies.item(0)

                ray = rays_sk.sketchCurves.sketchLines.addByTwoPoints(ZERO, point3d(INNERRADIUS * face.unitcenter))
                cplaneIn = cplanes.createInput()
                position = adsk.core.ValueInput.createByReal(1)
                cplaneIn.setByDistanceOnPath(ray, position)
                cplane = cplanes.add(cplaneIn)

                hole_sk = sketches.add(cplane)
                hole_sk.sketchCurves.sketchCircles.addByCenterRadius(ZERO, OUTERHOLE/2)
                hole_sk.sketchCurves.sketchCircles.addByCenterRadius(ZERO, INNERHOLE/2)

                innerholeprof = adsk.core.ObjectCollection.create()
                for profile in hole_sk.profiles:
                    innerholeprof.add(profile)

                hole_profiles = sorted(hole_sk.profiles, key=lambda profile: profile.areaProperties().area)
                outerholeprof = adsk.core.ObjectCollection.create()         
                outerholeprof.add(hole_profiles[1])

                extrudeIn = extrudes.createInput(outerholeprof, adsk.fusion.FeatureOperations.CutFeatureOperation)
                distance = adsk.core.ValueInput.createByReal(RADIUS)
                extrudeIn.setSymmetricExtent(distance, False)
                extrudeIn.participantBodies = facebodies
                extrudes.add(extrudeIn)

                extrudeIn = extrudes.createInput(innerholeprof, adsk.fusion.FeatureOperations.CutFeatureOperation)
                distance = adsk.core.ValueInput.createByReal(RADIUS)
                extrudeIn.setSymmetricExtent(distance, False)
                extrudeIn.participantBodies = corebodies
                extrudes.add(extrudeIn)

                progress4.update()

        app.log(h.getobstaclecount())


    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))


def point3d(vector):
    return adsk.core.Point3D.create(vector[0], vector[1], vector[2])

def point2d(x = None, y = None, r = None, theta = None, origin = None, rotation = math.pi/2):

    if x is None and y is None:
        x = r*math.cos(math.pi*theta/180)
        y = r*math.sin(math.pi*theta/180)

    X = y * math.cos(rotation) + x * math.sin(rotation)
    Y = y * math.sin(rotation) - x * math.cos(rotation)

    if origin:
        X = X + origin[0]
        Y = Y + origin[1]

    return adsk.core.Point3D.create(X, Y ,0)

class Progress:
    def __init__(self, ui, title, sequence):
        n = len(sequence)
        self.ui = ui
        self.dialog = ui.createProgressDialog()
        self.dialog.cancelButtonText = 'Give up'
        self.dialog.isCancelButtonShown = True
        self.dialog.show(title, "%p% complete (Step %v of %m)\nEstimated time remaining: ", 0, n)
        self.timer = Timer(n)

    def update(self):
        timestr = self.timer.report()
        self.dialog.progressValue = self.dialog.progressValue + 1
        self.dialog.message = f"%p% ({self.timer.percent}% adjusted)\nStep %v / %m\n{timestr}"

    def check(self):
        if self.dialog.wasCancelled:
            raise Exception("Script manually cancelled.")

    def finish(self):
        self.dialog.hide()


class Timer:
    def __init__(self, n):
        self.n = n

        self.t_0 = time.time()
        self.t_i = time.time()
        self.i = 0

        self.avgdt = 0
        self.idt = []

        (self.pt, self.ft, self.tt) = (0, math.inf, math.inf)
        
        self.r = []
        self.r_med = 1

        self.percent = 0

    def update(self):
        self.i = self.i + 1
        self.pt = time.time() - self.t_0
        self.avgdt = self.pt / self.i

        dt = time.time() - self.t_i
        k = self.n - self.i
        if self.idt == []:
            self.ft = k * dt
        else:
            r_i = dt / self.idt[-1]
            self.r.append(r_i)
            
            self.r_med = median(self.r)

            if r_i == 1:
                self.ft = self.ft - dt
            else:
                self.ft = self.avgdt * (1 - self.r_med**(k+1)) / (1 - self.r_med)


        self.tt = self.pt + self.ft
        self.percent = round(100 * self.pt / self.tt)
        self.idt.append(dt)
        self.t_i = time.time()

        return (self.pt, self.ft, self.tt)

    def timestring(self, t = None):
        if t == None:
            t = time.time() - self.t_0
        minutes, seconds = divmod(t, 60)
        hour, minutes = divmod(minutes, 60)
        return '%02d:%02d:%02d' % (hour, minutes, seconds)

    def report(self):
        self.update()

        strpt     = f"Elapsed time:                         {self.timestring(self.pt)}"
        strft     = f"Estimated time remaining:      {self.timestring(self.ft)}"
        strtt = f"Estimated total time:              {self.timestring(self.tt)}"
        strdt = f"Last step time: {self.idt[-1]:.3f} seconds"
        stravgdt = f"Average step time: {self.avgdt:.3f} seconds"

        return f"{strpt}\n{strft}\n{strtt}\n\n{strdt}\n{stravgdt}\nr = {self.r_med}"

    def total(self):
        return self.timestring(self.tt)
