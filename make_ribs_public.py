#coding:utf-8

__source__ =  'add算法建筑社'
__wechat__ = 'DA_Design'
__author__ = 'Forest Hell'


import rhinoscriptsyntax as rs
import System
import math
from itertools import combinations as combo
from functools import reduce as reduce
from Rhino import *
from Rhino.Input import *
from Rhino.DocObjects import *
from Rhino.Geometry import *
from scriptcontext import doc

def _caculate_thick(dir1, dir2, thickness):
    ''' determine cut width accroding to the angle between ribs'''
    theta = Vector3d.VectorAngle(dir1, dir2)
    theta2 = Vector3d.VectorAngle(-dir1, dir2)
    t = min(theta,theta2)
    return thickness*(1+math.cos(t)/math.sin(t))

def _cap_curve(c1, c2):
    line1 = LineCurve(c1.PointAtStart, c2.PointAtStart)
    line2 = LineCurve(c1.PointAtEnd, c2.PointAtEnd )
    return line1, line2

def _endp(curve):
    return curve.PointAtStart, curve.PointAtEnd

def cull_plane(planes):
    clean_planes = [planes[0]]
    for i in range(len(planes)):
        repeat = False
        for j in range(i+1,len(planes)):
            repeat = is_coplane(planes[i], planes[j])
        if not repeat:
            clean_planes.append(planes[i])
    return clean_planes

def is_coplane(pl1,pl2, epsilon = 0.1):
    return pl1.ZAxis.IsParallelTo(pl2.ZAxis, epsilon) and pl1.DistanceTo(pl2.Origin) <= epsilon

def offset(curve, distance, off_pl = None, cornor = CurveOffsetCornerStyle.Sharp):
    if not off_pl:
        off_pl = curve.TryGetPlane()[1]
    if hasattr(curve, '__iter__'):
        curve_list = []
        for i in curve:
            c = i.Offset(off_pl, distance,tol, cornor)
            if c: curve_list.append(c[0])
        return curve_list
    else:
        offset = curve.Offset(off_pl, distance,tol, cornor)[0]
        return offset

def align_curve(curve, plane, Orientation):
    if curve.IsClosed and curve.ClosedCurveOrientation(plane) != Orientation:
        curve.Reverse()

def render(objs):
    if hasattr(objs, '__iter__'):
        for obj in objs:
            render(obj)
    elif isinstance(objs, GeometryBase):
        id = doc.Objects.Add(objs)

def frame_curve(curve, step, epsilon = 1):
    '''perpendicular plane along curve'''
    domain = curve.Domain
    count = int(curve.GetLength() // step)
    frame = []
    for i in range(count):
        l = i * step
        t = domain[0] + curve.LengthParameter(l)[1]
        origin = curve.PointAt(t)
        normal = curve.TangentAt(t)
        this_pl = Plane(origin, normal)
        frame.append(this_pl)
    return cull_plane(frame)

def intsect(breps, thickness):
    combs = combo(breps, 2)
    for pair in combs:
        b1, b2 = pair
        res, curves, _ = Intersect.Intersection.BrepBrep(b1, b2, tol)
        if res:
            dir1 = b1.Surfaces[0].TryGetPlane()[1].Normal
            dir2 = b2.Surfaces[0].TryGetPlane()[1].Normal
            if dir1.IsParallelTo(dir2, 0.001) : continue
            width = _caculate_thick(dir1, dir2, thickness)
            cutlines[str(b1)].extend( make_notch(curves, width, b1, 1) )
            cutlines[str(b2)].extend( make_notch(curves, width, b2, -1) )

def contour(brep, curve, step, off_dist):
    '''make contour line along curve'''
    if off_dist >= 0:orient = CurveOrientation.Clockwise
    else:orient = CurveOrientation.CounterClockwise
    frames = frame_curve(curve, step)
    contour_curves = []
    offcurves = []
    Xcurve = []
    for plane in frames:
        if brep.ObjectType == ObjectType.Brep:
            res, Xcurve, _ = Intersect.Intersection.BrepPlane(brep, plane, tol)
        elif brep.ObjectType == ObjectType.Mesh:
            polylines = Intersect.Intersection.MeshPlane(brep, plane)
            if polylines:
                Xcurve = [pline.ToNurbsCurve() for pline in polylines]
        if Xcurve:
            Xcurve = Curve.JoinCurves(Xcurve)
            contour_curves.append(Xcurve)
            for c in Xcurve:
                align_curve(c, plane, orient)
            offcurves.append(offset(Xcurve,off_dist,plane))
    return contour_curves, offcurves

def make_ribs(brep, curve, step, off_dist, thickness):
    contour_curves, offcurves = contour(brep, curve, step, off_dist)
    edges = [ [contour_curves[i], offcurves[i]] for i in range(len(contour_curves)) ]
    # edges : [ [[curve1,curve2...],[offcurve1,offcurve2...]], ...]
    surf = []
    outlines = []
    for tup in edges:
        # 0: pack curve, pairs: [ [curve1, offcurve1], [curve12, offcurve2], ...]
        pairs = zip(tup[0], tup[1])
        pairs = [list(t) for t in pairs]
        unmerged_breps = []
        # 1: Create Planar Breps
        for pair in pairs:
            if not pair[0].IsClosed:
                pair.extend(_cap_curve(pair[0], pair[1]))
            # generate Non-hollow waffle if off_dist < 0
            if solid:unmerged_breps.extend( Brep.CreatePlanarBreps(pair[0]) )
            else:unmerged_breps.extend( Brep.CreatePlanarBreps(pair) )
        # 2: union interseted breps
        if not __author__ : raise
        merged = Brep.CreateBooleanUnion(unmerged_breps, tol) if len(unmerged_breps) > 1 else unmerged_breps
        surf.extend(merged)
        for brep in merged:
            inner_loops = Curve.JoinCurves(brep.DuplicateNakedEdgeCurves(False,True))
            outter_loops = Curve.JoinCurves(brep.DuplicateNakedEdgeCurves(True,False))
            # 3 : when brep is a loop , add extra cut line.if not ,inner_loops is empty
            outline = []
            extra_cuts = []
            for loop in inner_loops:
                pt = _endp(loop)[0]
                cpts = loop.ClosestPoints(outter_loops)
                
                extra_cuts.append(LineCurve(cpts[1], cpts[2]))
            outline = list(inner_loops + outter_loops)
            outline.extend(extra_cuts)
            outlines.append(outline)
    return surf, outlines

def get_next_len_pt(curve,pt,length):
    t = curve.ClosestPoint(pt)[1]
    i = Interval(t,curve.Domain[1])
    p = curve.LengthParameter(length, i)
    return curve.PointAt(p)
def make_notch(int_curves, width, brep, direct):
    ''' get notch rectangle '''
    dist = width / 2
    surface = brep.Faces[0]
    pl = surface.TryGetPlane()[1]
    notch = []
    for line in int_curves:
        off1 = offset(line, dist, pl)
        off2 = offset(line, -dist, pl)
        line1, line2 = _cap_curve(off1, off2)
        cut_rect = Curve.JoinCurves([off1,off2,line1,line2])[0]
        s, e = _endp(line)
        vec = Vector3d(e - s) / 2 if direct > 0 else  -Vector3d(e - s) / 2
        cut_rect.Translate(vec)
        notch.append(cut_rect)
    return notch


def main():
    #select a brep and 2 guide curves
    filter = ObjectType.Mesh | ObjectType.Brep
    result, obj = RhinoGet.GetOneObject('选择多重曲面或者网格', False, filter)
    brep = obj.Brep() or obj.Mesh()
    id1 = rs.GetCurveObject('第一条引导曲线')[0]
    curve1 = rs.coercecurve(id1)
    id2 = rs.GetCurveObject('第二条引导曲线')[0]
    while id2 == id1:
        id2 = rs.GetCurveObject('选择一条不同的曲线！')[0]
    curve2 = rs.coercecurve(id2)
    step1 = rs.GetInteger('沿第一条曲线的切割间距', number = 5)
    step2 = rs.GetInteger('沿第二条曲线的切割间距', number = 5)
    global solid
    solid = rs.GetBoolean('是否切穿实体？',('选择','否','是'),False)[0]
    offset_dist = rs.GetReal('生成条带的偏移宽度？（正值向外，负值向内,决定标注字体大小）', number = 2)
    thickness = rs.GetReal('材料的厚度', number = 0.3)

    #make ribs
    breps1, outline1 = make_ribs(brep, curve1, step1, offset_dist, thickness)
    breps2, outline2 = make_ribs(brep, curve2, step2, offset_dist, thickness)
    offset_dist = abs(offset_dist)
    #cutline stroage
    
    for s in breps1:
        cutlines[str(s)] = []
    for s in breps2:
        cutlines[str(s)] = []
    #find intersect and make notch line
    brep_surf = breps1 + breps2
    out_lines = outline1 + outline2
    intsect(brep_surf, thickness)
    #arrange to xy plane
    y1 = 0
    y2 = 0
    x_max = 0
    for i, b in enumerate(brep_surf):
        outline = out_lines[i]
        cut = cutlines[str(b)]
        for c in cut:
            count = 0
            for c2 in outline[:-1]:
                event = Intersect.Intersection.CurveCurve(c, c2, 0.2, tol)
                if event:
                    count += 1
            if count > 1:
                outline = outline[0:-1]
        #0: orient to worldXY plane
        surf = b.Surfaces[0]
        from_pl = surf.TryGetPlane()[1]
        to_pl = Plane.WorldXY
        xform = Transform.PlaneToPlane(from_pl, to_pl)
        group1 = Group(str(b))
        group1.Add(cut)
        group1.Add(outline)
        group = Group('new' + str(b))
        group.Add(group1.objects)
        b_in_breps1 = 0
        if b in breps1:
            text = 'x%d'%i
            b_in_breps1 = 1
        else:
            text = 'y%d'%(i-len(breps1))
        original_box = group1.bbox(plane = from_pl)
        group1.Add([TextDot(text,from_pl.Origin)])
        group1.Display()
        #arrange 1 : get boundingbox
        group.Transform(xform)
        bbox = group.bbox()
        bbox.Inflate(offset_dist, offset_dist, 0)
        xlen = (bbox.Max-bbox.Min).X
        ylen = (bbox.Max-bbox.Min).Y
        #arrange 2 : rotate according to x,y size
        if ylen > xlen:
            xlen, ylen = ylen, xlen
            rot = Transform.Rotation(math.pi/2, Vector3d.ZAxis,bbox.Center)
            group.Transform(rot)
        x_max = xlen if xlen > x_max else x_max
        #arrange 3 : align bbox cornor to worldXY origin
        cornor_plane = Plane(bbox.Min, Vector3d.ZAxis)
        align_xform = Transform.PlaneToPlane(cornor_plane, to_pl)
        group.Transform(align_xform)
        #arrange 4 : translate each group to proper position
        center = bbox.Center
        if b_in_breps1:
            move = Vector3d(0,y1,0)
            y1 += ylen 
        else:
            move = Vector3d(x_max ,y2,0)
            y2 += ylen 
        group.Transform(Transform.Translation(move))
        group.label(text,offset_dist)
        group.Display()
    render(brep_surf)
    doc.Objects.Hide(obj,True)
    doc.Views.Redraw()

class Group:
    def __init__(self, name):
        self.objects = []
        self.group = doc.Groups.Add(name)
        self.box = None
        self.annotation = []
    def bbox(self, plane = Plane.WorldXY):
        union = BoundingBox.Union
        boxes = [ geo.GetBoundingBox(plane) for geo in self.objects if isinstance(geo, GeometryBase)]
        lamb = lambda x, y: union(x, y)
        self.box = reduce(lamb, boxes)
        return self.box
    def Add(self, geo):
        self.objects.extend(geo)
    def Transform(self, xform):
        if self.box:
            self.box.Transform(xform)
        for obj in self.objects:
            obj.Transform(xform)
    def Display(self):
        for obj in self.objects:
            id = doc.Objects.Add(obj)
            doc.Groups.AddToGroup(self.group, id)
            doc.Groups.Add()
        for text in self.annotation:
            attr = ObjectAttributes()
            attr.LayerIndex = text_layer
            id = doc.Objects.Add(text,attr)
            doc.Groups.AddToGroup(self.group, id)
    def label(self, text, height,font = 'Arial'):
        center = self.box.Center
        t = TextEntity()
        t.FontIndex = 0
        t.Text = text
        t.TextHeight = height
        if not __source__ : raise
        X = Vector3d(self.box.Corner(False,True,True) - self.box.Corner(True,True,True))
        Y = Vector3d(self.box.Corner(True,False,True) - self.box.Corner(True,True,True))
        t.Plane = Plane(center,X, Y)
        text_curve = Curve.JoinCurves(t.Explode())
        self.annotation.extend(text_curve)
    def Add_id(self,id):
        doc.Groups.AddToGroup(self.group, id)

solid = False
text_layer = doc.Layers.Add('text',System.Drawing.Color.Red)
if text_layer == -1:
    text_layer = doc.Layers.Find('text',True)
tol = doc.PageAbsoluteTolerance
cutlines = {}
main()
