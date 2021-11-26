import math
import json
import jsonschema

# Hetool Library
# Half-Edge Based Data Structure for Two-Dimensional Solid Modeling
# Main developer: Danilo Silva Bomfim (dsbomfim2@hotmail.com)
# Contributors:
#              - André M. B. Pereira
#              - Luiz F. Bez
#              - Luiz F. Martha
#              - Pedro C. F. Lopes
#              - Rodrigo L. Soares

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# --------------------------------- GEOMETRY ---------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


class Point():

    def __init__(self, _x=None, _y=None):
        self.x = _x
        self.y = _y
        self.selected = False
        self.vertex = None
        self.attributes = []

    def setX(self, _x):
        self.x = _x

    def setY(self, _y):
        self.y = _y

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def setCoords(self, _x, _y):
        self.x = _x
        self.y = _y

    def setSelected(self, _select):
        self.selected = _select

    def isSelected(self):
        return self.selected

    # Equality test with tolerance (Manhattan distance)
    @staticmethod
    def equal(p1, p2, tol):
        return abs(p1.x - p2.x) < tol.x and abs(p1.y-p2.y) < tol.y

    # Equality test without tolerance
    def __eq__(p1, p2):
        return (p1.x == p2.x) and (p1.y == p2.y)

    # operator <
    def __lt__(p1, p2):
        if p1.x == p2.x:
            return p1.y < p2.y
        else:
            return p1.x < p2.x

    # # operator >
    def __gt__(p1, p2):
        if p1.x == p2.x:
            return p1.y > p2.y
        else:
            return p1.x > p2.x

    # Inequality test without tolerance
    def __ne__(p1, p2):
        return not (p1 == p2)

    # Addition +
    def __add__(p1, p2):
        return Point(p1.x+p2.x, p1.y+p2.y)

    # Addition +=
    def __iadd__(p1, p2):
        p1 = p1+p2
        return p1

    # Subtraction -
    def __sub__(p1, p2):
        return Point(p1.x-p2.x, p1.y-p2.y)

    # Subtraction -=
    def __isub__(p1, p2):
        p1 = p1-p2
        return p1

    # Scalar multiplication
    def __mul__(p, s):
        return Point(p.x*s, p.y*s)

    # Division by scalar
    def __truediv__(p, s):
        if s == 0:
            return Point(0.0, 0.0)
        return Point(p.x/s, p.y/s)

    # Euclidian distance
    @staticmethod
    def euclidiandistance(p1, p2):
        return math.sqrt((p1.x-p2.x)*(p1.x-p2.x) +
                         (p1.y-p2.y)*(p1.y-p2.y))

    # Manhattan distance
    @staticmethod
    def manhattandistance(p1, p2):
        return abs(p1.x-p2.x) + abs(p1.y-p2.y)

    # Square of size (vector)
    @staticmethod
    def sizesquare(p):
        return (p.x*p.x + p.y*p.y)

    # Size (vector)
    @staticmethod
    def size(p):
        return math.sqrt(Point.sizesquare(p))

    # Dot product (vector)
    @staticmethod
    def dotprod(p1, p2):
        return p1.x*p2.x + p1.y*p2.y

    # Cross product (vector)
    @staticmethod
    def crossprod(p1, p2):
        return p1.x*p2.y - p2.x*p1.y

    # Normalize (vector)
    @staticmethod
    def normalize(p):
        norm = Point.size(p)
        if norm == 0:
            return Point(0.0, 0.0)
        return Point(p.x/norm, p.y/norm)

    # Twice the signed area of triangle p1-p2-p3
    @staticmethod
    def area2d(p1, p2, p3):
        ptA = p2 - p1
        ptB = p3 - p1
        return Point.crossprod(ptA, ptB)


class Segment:
    selected = False
    PARAM_TOL = 1e-7
    nsudv = None

    def setNumberOfSubdivisions(self, _number):
        self.nsudv = _number

    def getNumberOfSubdivisions(self):
        return self.nsudv

    def setSelected(self, _select):
        self.selected = _select

    def isSelected(self):
        return self.selected


class Line(Segment):

    def __init__(self, _pt1=None, _pt2=None):
        self.pt1 = _pt1
        self.pt2 = _pt2
        self.nPts = 0
        self.edge = None
        self.attributes = []

        if _pt1 is not None:
            self.nPts += 1

        if _pt2 is not None:
            self.nPts += 1

    def addPoint(self, _x, _y):
        if self.nPts == 0:
            self.pt1 = Point(_x, _y)
            self.nPts += 1
        elif self.nPts == 1:
            self.pt2 = Point(_x, _y)
            self.nPts += 1

    def getNumberOfPoints(self):
        return self.nPts

    def getPoint(self, _t):
        vx = self.pt2.getX() - self.pt1.getX()
        vy = self.pt2.getY() - self.pt1.getY()
        if _t < 0:
            xOn = self.pt1.getX()
            yOn = self.pt1.getY()
        elif _t > 1:
            xOn = self.pt2.getX()
            yOn = self.pt2.getY()
        else:
            xOn = self.pt1.getX() + _t * vx
            yOn = self.pt1.getY() + _t * vy
        return Point(xOn, yOn)

    def isPossible(self):
        if self.nPts < 2:
            return False
        return True

    def getPoints(self):
        tempPts = []
        if self.nPts == 1:
            tempPts.append(self.pt1)
            return tempPts

        tempPts.append(self.pt1)
        tempPts.append(self.pt2)
        return tempPts

    def getPointsToDraw(self):
        tempPts = []
        tempPts.append(self.pt1)
        tempPts.append(self.pt2)
        return tempPts

    def getPointsToDrawPt(self, _pt):
        tempPts = []
        tempPts.append(self.pt1)
        if self.nPts == 2:
            tempPts.append(self.pt2)
        elif self.nPts == 1:
            tempPts.append(_pt)
        return tempPts

    def setInitPoint(self, _pt):
        self.pt1 = _pt

    def setEndPoint(self, _pt):
        self.pt2 = _pt

    def closestPoint(self, _x, _y):
        vx = self.pt2.getX() - self.pt1.getX()
        vy = self.pt2.getY() - self.pt1.getY()
        t = (vx*(_x - self.pt1.getX()) + vy *
             (_y - self.pt1.getY())) / (vx*vx + vy*vy)

        if t < 0.0:
            xOn = self.pt1.getX()
            yOn = self.pt1.getY()
        elif t > 1.0:
            xOn = self.pt2.getX()
            yOn = self.pt2.getY()
        else:
            xOn = self.pt1.getX() + t * vx
            yOn = self.pt1.getY() + t * vy

        dist = math.sqrt((xOn - _x)*(xOn - _x)+(yOn - _y)*(yOn - _y))
        return xOn, yOn, dist

    def getBoundBox(self):

        xmax = max(self.pt1.getX(), self.pt2.getX())
        xmin = min(self.pt1.getX(), self.pt2.getX())
        ymax = max(self.pt1.getY(), self.pt2.getY())
        ymin = min(self.pt1.getY(), self.pt2.getY())
        return xmin, xmax, ymin, ymax

    def getType(self):
        return 'LINE'

    def isUnlimited(self):
        return False

    def getXinit(self):
        return self.pt1.getX()

    def getYinit(self):
        return self.pt1.getY()

    def getXend(self):
        return self.pt2.getX()

    def getYend(self):
        return self.pt2.getY()

    def boundIntegral(self):
        return (self.pt1.getX()*self.pt2.getY() - self.pt2.getX()*self.pt1.getY())*0.5

    def length(self, _t0, _t1):
        p1 = self.getPoint(_t0)
        p2 = self.getPoint(_t1)
        len = math.sqrt((p2.getX() - p1.getX()) *
                        (p2.getX() - p1.getX()) +
                        (p2.getY() - p1.getY()) *
                        (p2.getY() - p1.getY()))

        return len

    def tangent(self, _t):
        pts = self.getPoints()
        tan = pts[1] - pts[0]
        tan = Point.normalize(tan)

        return tan

    def curvature(self, _t):
        return 0.0

    def selfIntersect(self):
        return False, None, None

    def clone(self):
        myClone = Line(self.pt1, self.pt2)
        return myClone

    def splitSegment(self, _t, _pt):
        if _t <= Segment.PARAM_TOL:
            _segment1 = None
            _segment2 = self
            return _segment1, _segment2

        if 1.0 - _t <= Segment.PARAM_TOL:
            _segment1 = self
            _segment2 = None
            return _segment1, _segment2

        _segment1 = Line(self.pt1, _pt)
        _segment2 = Line(_pt, self.pt2)

        return _segment1, _segment2

    def split(self, _params, _pts):
        seg2 = self.clone()
        segments = []

        for i in range(0, len(_params)):
            seg1, seg2 = seg2.splitSegment(_params[i], _pts[i])
            segments.append(seg1)

            # update the remaining parameters
            for j in range(i+1, len(_params)):
                _params[j] = (_params[j] - _params[i])/(1-_params[i])

        segments.append(seg2)

        return segments

    def intersectPoint(self, _pt, _tol):
        pts = self.getPoints()
        dist, pi, t = CompGeom.getClosestPointSegment(pts[0], pts[1], _pt)

        if dist <= _tol:
            return True, t, pi

        return False, t, pi

    def intersectSegment(self, _segment):

        poly = _segment.getPoints()
        if _segment.getType() == 'LINE':
            return CompGeom.computeLineIntersection(self.pt1, self.pt2, poly[0], poly[1])
        elif _segment.getType() == 'POLYLINE':
            pts = self.getPoints()
            return CompGeom.computePolyPolyIntersection(pts, poly)

    def isEqual(self, _segment, _tol):

        if _segment.getType() == 'LINE':

            pts1 = self.getPoints()
            pts2 = _segment.getPoints()
            tol = Point(_tol, _tol)

            if Point.equal(pts1[0], pts2[0], tol):
                if Point.equal(pts1[1], pts2[1], tol):
                    return True
            elif Point.equal(pts1[1], pts2[0], tol):
                if Point.equal(pts1[0], pts2[1], tol):
                    return True
            return False

        else:
            return _segment.isEqual(self, _tol)

    def ray(self, _pt):
        x = _pt.getX()
        y = _pt.getY()

        if self.pt1.getY() == self.pt2.getY():  # discard horizontal line
            return 0.0

        if self.pt1.getY() > y and self.pt2.getY() > y:  # discard line above ray
            return 0.0

        if self.pt1.getY() < y and self.pt2.getY() < y:  # Discard line below ray
            return 0.0

        if self.pt1.getX() < x and self.pt2.getX() < x:  # discard line to the left of point
            return 0.0

        if self.pt1.getY() == y:  # ray passes at first line point
            if self.pt1.getX() > x and self.pt2.getY() > y:
                # Intersects if first point is to the right of given point
                # and second point is above.
                return 1
        else:
            if self.pt2.getY() == y:  # ray passes at second point
                if self.pt2.getX() > x and self.pt1.getY() > y:
                    # Intersects if first point is to the right of given point
                    # and second point is above.
                    return 1
            else:  # ray passes with first and second points
                if self.pt1.getX() > x and self.pt2.getX() > x:
                    # Intersects if first point is to the right of given point
                    # and second point is above.
                    return 1
                else:
                    # Compute x coordinate of intersection of ray with line segment
                    dx = self.pt1.getX() - self.pt2.getX()
                    xc = self.pt1.getX()

                    if dx != 0:
                        xc += (y - self.pt1.getY())*dx / \
                            (self.pt1.getY() - self.pt2.getY())

                    if xc > x:
                        # Intersects if first point is to the right of given point
                        # and second point is above.
                        return 1

        return 0.0


class Polyline(Segment):
    def __init__(self, _pts=None):
        self.pts = _pts
        if self.pts is None:
            self.pts = []
        self.nPts = 0
        self.edge = None
        self.attributes = []

    def addPoint(self, _x, _y):
        self.pts.append(Point(_x, _y))
        self.nPts += 1

    def getNumberOfPoints(self):
        return self.nPts

    def getPoint(self, _t):

        if _t <= 0:
            return Point(self.pts[0].getX(), self.pts[0].getY())

        if _t >= 1.0:
            return Point(self.pts[-1].getX(), self.pts[-1].getY())

        length = self.length(0, 1)
        s = _t*length
        loc_t = 1.0
        prev_id = 0
        next_id = 0
        length = 0

        for i in range(1, len(self.pts)):
            prev_id = i - 1
            next_id = i
            dist = math.sqrt((self.pts[i].getX() - self.pts[i - 1].getX()) *
                             (self.pts[i].getX() - self.pts[i - 1].getX()) +
                             (self.pts[i].getY() - self.pts[i - 1].getY()) *
                             (self.pts[i].getY() - self.pts[i - 1].getY()))

            if (length + dist) >= s:
                loc_t = (s - length) / dist
                break

            length += dist

        x = self.pts[prev_id].getX() + loc_t * \
            (self.pts[next_id].getX() - self.pts[prev_id].getX())
        y = self.pts[prev_id].getY() + loc_t * \
            (self.pts[next_id].getY() - self.pts[prev_id].getY())

        return Point(x, y)

    def isPossible(self):
        if self.nPts < 2:
            return False

        return True

    def getPoints(self):
        return self.pts

    def getPointsToDraw(self):
        return self.pts

    def getPointsToDrawPt(self, _pt):
        tempPts = []
        for i in range(0, self.nPts):
            tempPts.append(self.pts[i])

        tempPts.append(_pt)
        return tempPts

    def setInitPoint(self, _pt):
        self.pts[0] = _pt

    def setEndPoint(self, _pt):
        self.pts[-1] = _pt

    def closestPoint(self, _x, _y):

        aux = Line(self.pts[0], self.pts[1])
        x, y, d = aux.closestPoint(_x, _y)
        xOn = x
        yOn = y

        dmin = d

        for i in range(2, len(self.pts)):
            aux = Line(self.pts[i - 1], self.pts[i])
            x, y, d = aux.closestPoint(_x, _y)

            if d < dmin:
                xOn = x
                yOn = y
                dmin = d

        return xOn, yOn, dmin

    def getBoundBox(self):
        x = []
        y = []
        for point in self.pts:
            x.append(point.getX())
            y.append(point.getY())

        xmin = min(x)
        xmax = max(x)
        ymin = min(y)
        ymax = max(y)

        return xmin, xmax, ymin, ymax

    def getType(self):
        return 'POLYLINE'

    def isUnlimited(self):
        return True

    def getXinit(self):
        return self.pts[0].getX()

    def getYinit(self):
        return self.pts[0].getY()

    def getXend(self):
        return self.pts[-1].getX()

    def getYend(self):
        return self.pts[-1].getY()

    def boundIntegral(self):
        area = 0

        for i in range(0, len(self.pts)-1):
            pt1 = self.pts[i]
            pt2 = self.pts[i+1]
            area += (pt1.getX())*(pt2.getY()) - (pt2.getX())*(pt1.getY())

        return area*0.5

    def curvature(self, _t):
        return 0.0

    def tangent(self, _t):

        if _t <= 0.0:
            tan = self.pts[1]-self.pts[0]
            tan = Point.normalize(tan)
            return tan

        if _t >= 1.0:
            tan = self.pts[-1]-self.pts[-2]
            tan = Point.normalize(tan)
            return tan

        length = 0.0

        for j in range(1, len(self.pts)):
            length += math.sqrt((self.pts[j].getX()-self.pts[j-1].getX())
                                * (self.pts[j].getX()-self.pts[j-1].getX()) +
                                (self.pts[j].getY()-self.pts[j-1].getY()) *
                                (self.pts[j].getY()-self.pts[j-1].getY()))

        s = _t*length
        prev_id = 0
        next_id = 0

        for i in range(1, len(self.pts)):
            prev_id = i - 1
            next_id = i

            d = math.sqrt((self.pts[i].getX()-self.pts[i-1].getX())
                          * (self.pts[i].getX()-self.pts[i-1].getX()) +
                          (self.pts[i].getY()-self.pts[i-1].getY()) *
                          (self.pts[i].getY()-self.pts[i-1].getY()))

            if length + d >= s:
                break

            length += d

        tan = self.pts[next_id] - self.pts[prev_id]
        tan = Point.normalize(tan)

        return tan

    def selfIntersect(self):
        flag, pts, params = CompGeom.splitSelfIntersected(self.getPoints())
        return flag, pts, params

    def clone(self):
        myClone = Polyline(self.pts)
        return myClone

    def length(self, _t0, _t1):
        L = 0.0
        pts = self.getPoints()
        for i in range(1, len(pts)):
            L += math.sqrt((pts[i].getX()-pts[i-1].getX())*(pts[i].getX()-pts[i-1].getX()) + (
                pts[i].getY()-pts[i-1].getY()) * (pts[i].getY()-pts[i-1].getY()))

        return L*(_t1-_t0)

    def splitSegment(self, _t, _pt):
        if _t <= Segment.PARAM_TOL:
            _segment1 = None
            _segment2 = self
            return _segment1, _segment2

        if 1.0 - _t <= Segment.PARAM_TOL:
            _segment1 = self
            _segment2 = None
            return _segment1, _segment2

        L = self.length(0, 1)
        s = _t*L
        loc_t = 1.0
        prev_id = 0
        next_id = 0
        L = 0.0
        pts = self.getPoints()
        for j in range(1, len(pts)):
            prev_id = j - 1
            next_id = j
            d = math.sqrt((pts[j].getX()-pts[j-1].getX())*(pts[j].getX()-pts[j-1].getX()) + (
                pts[j].getY()-pts[j-1].getY()) * (pts[j].getY()-pts[j-1].getY()))

            if (L+d) >= s:
                loc_t = (s-L)/d
                break
            L += d

        segment1_pts = []
        segment2_pts = []

        for i in range(0, prev_id):
            segment1_pts.append(pts[i])

        # check whether the split point is one of the points of the polyline itself
        if loc_t > Segment.PARAM_TOL:
            segment1_pts.append(pts[prev_id])

        segment1_pts.append(_pt)
        segment2_pts.append(_pt)

        if 1.0 - loc_t > Segment.PARAM_TOL:
            segment2_pts.append(pts[next_id])

        for j in range(next_id+1, len(pts)):
            segment2_pts.append(pts[j])

        _segment1 = Polyline(segment1_pts)

        _segment2 = Polyline(segment2_pts)

        return _segment1, _segment2

    def split(self, _params, _pts):
        seg2 = self.clone()
        segments = []

        for i in range(0, len(_params)):
            seg1, seg2 = seg2.splitSegment(_params[i], _pts[i])
            segments.append(seg1)

            # update the remaining parameters
            for j in range(i+1, len(_params)):
                _params[j] = (_params[j] - _params[i])/(1-_params[i])

        segments.append(seg2)

        return segments

    def ray(self, _pt):
        x = _pt.getX()
        y = _pt.getY()
        n = len(self.pts)
        ni = 0

        for i in range(0, n-1):
            pt1 = self.pts[i]
            pt2 = self.pts[i+1]

            if pt1.getY() == pt2.getY():  # discard horizontal line
                continue

            if pt1.getY() > y and pt2.getY() > y:  # discard line above ray
                continue

            if pt1.getY() < y and pt2.getY() < y:  # Discard line below ray
                continue

            if pt1.getX() < x and pt2.getX() < x:  # discard line to the left of point
                continue

            if pt1.getY() == y:  # ray passes at first line point
                if pt1.getX() > x and pt2.getY() > y:
                    # Count intersection if first point is to the right of given point
                    # and second point is above.
                    ni += 1
            else:
                if pt2.getY() == y:  # ray passes at second point
                    if pt2.getX() > x and pt1.getY() > y:
                        # Count intersection if first point is to the right of given point
                        # and second point is above.
                        ni += 1
                else:  # ray passes with first and second points
                    if pt1.getX() > x and pt2.getX() > x:
                        # Count intersection if first point is to the right of given point
                        # and second point is above.
                        ni += 1
                    else:
                        # Compute x coordinate of intersection of ray with line segment
                        dx = pt1.getX() - pt2.getX()
                        xc = pt1.getX()

                        if dx != 0:
                            xc += (y - pt1.getY())*dx / \
                                (pt1.getY() - pt2.getY())

                        if xc > x:
                            # Count intersection if first point is to the right of given point
                            # and second point is above.
                            ni += 1

        return ni

    def isEqual(self, _segment, _tol):
        tol = Point(_tol, _tol)
        if _segment.getType() == 'LINE':
            if len(self.pts) == 2:
                ptsLine = _segment.getPoints()
                if Point.equal(ptsLine[0], self.pts[0], tol):
                    if Point.equal(ptsLine[1], self.pts[1], tol):
                        return True

                elif Point.equal(ptsLine[0], self.pts[1], tol):
                    if Point.equal(ptsLine[1], self.pts[0], tol):
                        return True

                return False

            else:
                return False

        elif _segment.getType() == 'POLYLINE':
            thatPts = _segment.getPoints()

            if len(self.pts) != len(thatPts):
                return False

            if Point.equal(self.pts[0], thatPts[0], tol):
                for i in range(1, len(self.pts)):
                    if not Point.equal(self.pts[i], thatPts[i], tol):
                        return False

                return True

            else:

                for i in range(0, len(self.pts)):
                    if not Point.equal(self.pts[-1-i], thatPts[i], tol):
                        return False

                return True

        else:
            return _segment.isEqual(self, _tol)

    def intersectPoint(self, _pt, _tol):

        partialLength = 0
        totalLength = self.length(0, 1)
        interStatus = False
        param = None

        if Point.euclidiandistance(_pt, self.pts[0]) <= _tol:
            param = 0

        for i in range(1, len(self.pts)):
            p1 = Point(self.pts[i-1].getX(), self.pts[i-1].getY())
            p2 = Point(self.pts[i].getX(), self.pts[i].getY())

            dist, pi, t = CompGeom.getClosestPointSegment(p1, p2, _pt)
            length = math.sqrt(((self.pts[i].getX() - self.pts[i-1].getX()) *
                                (self.pts[i].getX() - self.pts[i-1].getX())) +
                               ((self.pts[i].getY() - self.pts[i-1].getY()) *
                                (self.pts[i].getY() - self.pts[i-1].getY())))

            # skip init intersections at each segment (no repeated intersections)
            if dist <= _tol and t*length > _tol:
                param = ((partialLength + t*length) / totalLength)
                interStatus = True
                break

            partialLength += length

        return interStatus, param, pi

    def intersectSegment(self, _segment):

        status = []
        param1 = []
        param2 = []
        pts = []

        if _segment.getType() == 'LINE' or _segment.getType() == 'POLYLINE':
            poly = _segment.getPoints()
            return CompGeom.computePolyPolyIntersection(self.pts, poly)
        else:
            # for each segment, create a line and intersect each line with the given segment
            totalLength = 0.0
            segLength = 0.0

            for i in range(0, len(self.pts)-1):
                segment = Line(self.pts[i], self.pts[i + 1])
                segLength = segment.length(0, 1)
                segPts, segmentParams, segParams = _segment.intersectsegment(
                    segment)

                for i in range(0, len(segPts)):
                    pts.append(segPts[i])
                    param1.append((segParams[i]*segLength) + totalLength)
                    param2.append(segmentParams[i])

                status.append(len(param1) > 0)
                totalLength += segLength
                segParams = []
                segmentParams = []
                segPts = []

        return status, pts, param1, param2


class Patch:

    def __init__(self):
        self.pts = []  # boundary points
        self.segments = []  # vector of boundary segments
        # orientations of segments with respect to counter-clockwise region boundary
        self.segmentOrients = []
        self.mesh = None
        self.selected = False
        self.holes = []  # vector of region holes
        self.holesOrients = []
        self.internalSegments = []
        self.internalSegmentsOrients = []
        self.isDeleted = False
        self.face = None
        self.attributes = []

    def __del__(self):
        if self.mesh:
            del self.mesh

    def getPoints(self):
        return self.pts

    def getSegments(self):
        return self.segments

    def getSegmentOrients(self):
        return self.segmentOrients

    def setSelected(self, _select):
        self.selected = _select

    def isSelected(self):
        return self.selected

    def setMesh(self, _mesh):
        self.mesh = _mesh

    def getMesh(self):
        return self.mesh

    def getBoundBox(self):

        if len(self.pts) == 0:
            return

        xmin = self.pts[0].getX()
        ymin = self.pts[0].getY()
        xmax = self.pts[0].getX()
        ymax = self.pts[0].getY()

        if len(self.pts) == 1:
            return

        for j in range(1, len(self.pts)):
            xmin = min(xmin, self.pts[j].getX())
            xmax = max(xmax, self.pts[j].getX())
            ymin = min(ymin, self.pts[j].getY())
            ymax = max(ymax, self.pts[j].getY())

        return xmin, xmax, ymin, ymax

    def setBoundary(self, _boundarysegments, _isOriented):
        self.segments = _boundarysegments.copy()
        self.segmentOrients = _isOriented.copy()
        self.pts = self.boundaryPolygon()

    def setHoles(self, _holessegments, _isOriented):
        self.holes = _holessegments
        self.holesOrients = _isOriented

    def setInternalSegments(self, _internalSegments, _isOriented):
        self.internalSegments = _internalSegments
        self.internalSegmentsOrients = _isOriented

    def isPointInside(self, _pt):
        numIntersec = 0
        for i in range(0, len(self.segments)):
            numIntersec += self.segments[i].ray(_pt)

        if numIntersec % 2 != 0:
            for i in range(0, len(self.holes)):
                numIntersec = 0
                for j in range(0, len(self.holes[i])):
                    numIntersec += self.holes[i][j].ray(_pt)

                if numIntersec % 2 != 0:
                    return False

            return True

        else:
            return False

    def boundaryPolygon(self):
        polygon = []
        for i in range(0, len(self.segments)):
            segmentPol = self.segments[i].getPoints()
            if self.segmentOrients[i]:
                for j in range(0, len(segmentPol)-1):
                    polygon.append(segmentPol[j])
            else:
                for j in range(len(segmentPol)-1, 0, -1):
                    polygon.append(segmentPol[j])

        return polygon

    def boundaryHole(self):
        polygons = []

        for i in range(0, len(self.holes)):
            polygon = []
            for j in range(0, len(self.holes[i])):
                segmentpol = self.holes[i][j].getPoints()
                if self.holesOrients[i][j]:
                    for m in range(0, len(segmentpol)-1):
                        polygon.append(segmentpol[m])
                else:
                    for m in range(len(segmentpol)-1, 0, -1):
                        polygon.append(segmentpol[m])

            polygon.reverse()
            polygons.append(polygon)

        return polygons

    def boundaryInternalSegments(self):
        polygons = []
        for i in range(0, len(self.internalSegments)):
            polygon = []
            for j in range(0, len(self.internalSegments[i])):
                segmentpol = self.internalSegments[i][j].getPoints()
                if self.internalSegmentsOrients[i][j]:
                    for m in range(0, len(segmentpol)-1):
                        polygon.append(segmentpol[m])
                else:
                    for m in range(len(segmentpol)-1, 0, -1):
                        polygon.append(segmentpol[m])

            polygon.reverse()
            polygons.append(polygon)

        return polygons

    def Area(self):
        Area = 0
        pts = self.pts
        triangs = Tesselation.triangleParing(pts)
        for j in range(0, len(triangs)):
            a = Point(pts[triangs[j][0]].getX(),
                      pts[triangs[j][0]].getY())
            b = Point(pts[triangs[j][1]].getX(),
                      pts[triangs[j][1]].getY())
            c = Point(pts[triangs[j][2]].getX(),
                      pts[triangs[j][2]].getY())

            Area += (a.getX()*b.getY() - a.getY()*b.getX()
                     + a.getY()*c.getX() - a.getX()*c.getY()
                     + b.getX()*c.getY() - c.getX()*b.getY()) / 2.0

        internalFaces = self.face.internalFaces()
        if len(internalFaces) > 0:
            for face in internalFaces:
                Area -= face.patch.Area()
                adjacentFaces = face.adjacentFaces()
                for adjface in adjacentFaces:
                    if adjface not in internalFaces and adjface != self.face:
                        pts = adjface.patch.getPoints()
                        if CompGeom.isPointInPolygon(self.pts, pts[0]):
                            Area -= adjface.patch.Area()

        return Area


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ---------------------------- ATTRIBUTE MANAGER -----------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


class AttribManager:

    def __init__(self):
        self.prototypes = []
        self.attributes = []

        with open("attribprototype.json", 'r') as file:
            input = json.load(file)

        with open("attrib_schema.json", 'r') as file:
            schemas = json.load(file)

        self.prototypes = input['prototypes']

        invalid_prototypes = []
        # valid attribute prototypes
        for prototype in self.prototypes:
            check = AttribManager.validate_attribute(prototype, schemas)
            if not check:
                invalid_prototypes.append(prototype)

        # remove invalid attribute prototypes
        for att in invalid_prototypes:
            self.prototypes.remove(att)

    def getPrototypes(self):
        return self.prototypes

    def getAttributes(self):
        return self.attributes

    def getAttributeByName(self, _name):
        for attribute in self.attributes:
            if attribute['name'] == _name:
                return attribute

    def removeAttribute(self, _attribute):
        self.attributes.remove(_attribute)

    def getPrototypeByType(self, _type):
        for prototype in self.prototypes:
            if prototype['type'] == _type:
                return prototype

    def createAttributeFromPrototype(self, _type, _name):

        # checks if an attribute with that name already exists
        for atributte in self.attributes:
            if atributte['name'] == _name:
                return False

        # get the prototype to create the name of the new attribute
        for prototype in self.prototypes:
            if prototype['type'] == _type:
                prototype_target = prototype.copy()
                prototype_target['properties'] = prototype['properties'].copy()
                self.attributes.append(prototype_target)
                prototype_target['name'] = _name
                return True

    def setAttributeValues(self, _name, _values):

        attribute = self.getAttributeByName(_name)
        attValues = attribute['properties']

        index = 0
        for key in attValues:
            attValues[key] = _values[index]
            index += 1

    @staticmethod
    def validate_attribute(_attribute, _schemas):

        check = AttribManager.validate_schema(
            _attribute, _schemas['att_schema'])
        if not check:
            return False

        values_types = []
        properties = _attribute['properties']
        for key in properties:
            valueType = type(properties[key])
            if valueType == int:
                values_types.append("int")
            elif valueType == float:
                values_types.append("float")
            elif valueType == bool:
                values_types.append("bool")
            elif valueType == str:
                values_types.append("string")
            elif valueType == list:
                values_types.append("color")
                color_dict = properties[key]
                check = AttribManager.validate_schema(
                    color_dict, _schemas['color_schema'])
                if not check:
                    return False
            elif valueType == dict:
                values_types.append("options")
                options_dict = properties[key]
                check = AttribManager.validate_schema(
                    options_dict, _schemas['options_schema'])
                if not check:
                    return False

                if options_dict['index'] < 0 or options_dict['index'] > (len(options_dict['list'])-1):
                    print("Given Attribute is InValid: index out of range")
                    print(_attribute)
                    return False

        # creates a new key
        _attribute['properties_type'] = values_types.copy()

        return True

    @staticmethod
    def validate_schema(_objetc, _schema):
        try:
            jsonschema.validate(instance=_objetc, schema=_schema)
            return True
        except jsonschema.exceptions.ValidationError as err:
            print("Given Attribute is InValid")
            print(_objetc)
            print(err)
            return False


class AttribSymbols:

    @staticmethod
    def getSymbol(_attribute, _scale, _pt=None, _seg=None, _patch=None):

        lines = []
        triangles = []
        squares = []
        points = []
        circles = []
        time = "before"  # determines whether the attribute will be drawn before or after

        if _attribute['symbol'] == 'Support':
            time = "before"
            if _pt is not None:
                lines, triangles, squares, circles = AttribSymbols.supportPoint(
                    _attribute, _pt, _scale)
            else:
                lines, triangles, squares, circles = AttribSymbols.supportSegment(
                    _attribute, _seg, _scale)

        elif _attribute['symbol'] == 'Arrow':

            if _attribute['type'] == "Concentrated Load":
                time = "after"
                if _pt is not None:
                    lines, triangles, circles = AttribSymbols.arrowPointCL(
                        _attribute, _pt, _scale)

            elif _attribute['type'] == "Uniform Load":
                time = "after"
                if _seg is not None:
                    lines, triangles = AttribSymbols.arrowSegmentUL(
                        _attribute, _seg, _scale)

        elif _attribute['symbol'] == 'Nsbdvs':
            time = "after"
            if _seg is not None:
                points = AttribSymbols.Nsbdvs(_attribute, _seg)

        # get the colors
        colors = []
        index = 0
        for att_type in _attribute['properties_type']:
            if att_type == "color":
                colors.append(list(_attribute['properties'].values())[
                    index].copy())
            index += 1

        symbol = {
            "lines": lines,
            "triangles": triangles,
            "squares": squares,
            "circles": circles,
            "colors": colors,
            "points": points,
            "time": time
        }

        return symbol

    @ staticmethod
    def rotateCoord(_pt, _ang):
        pt = Point(_pt.getX(), _pt.getY())
        x = (pt.x*(math.cos(_ang))) + (pt.y*(math.sin(_ang)))
        y = (pt.y*(math.cos(_ang))) - (pt.x*(math.sin(_ang)))
        return Point(x, y)

    @ staticmethod
    def getAngWithXDirec(_v2):
        v1 = Point(1, 0)
        ang = math.acos(Point.dotprod(v1, _v2)/((Point.size(v1)) *
                                                (Point.size(_v2))))
        ang = ang*180/CompGeom.PI

        return ang

    @ staticmethod
    def triangleSymbol(_pt, _scale, _ang):

        _ang = _ang*CompGeom.PI/180
        x = AttribSymbols.rotateCoord(Point(1*_scale, 0), _ang)
        y = AttribSymbols.rotateCoord(Point(0, 1*_scale), _ang)

        pt_a = _pt - x*0.75
        pt_b = pt_a + (y/2)
        pt_c = pt_a - (y/2)

        return [_pt, pt_b, pt_c]

    def squareSymbol(_pt, _scale, _ang):

        _ang = _ang*CompGeom.PI/180
        x = AttribSymbols.rotateCoord(Point(1*_scale, 0), _ang)
        y = AttribSymbols.rotateCoord(Point(0, 1*_scale), _ang)

        pt_a = _pt - (x/4) + (y/4)
        pt_b = _pt + (x/4) + (y/4)
        pt_c = _pt + (x/4) - (y/4)
        pt_d = _pt - (x/4) - (y/4)

        return [pt_a, pt_b, pt_c, pt_d]

    def circleSymbol(_pt, _r):
        x = _pt.getX()
        y = _pt.getY()
        num = 30
        circ_points = []

        for i in range(0, num):
            theta = 2*CompGeom.PI*i/num
            pt = Point(x + _r*math.cos(theta), y + _r*math.sin(theta))
            circ_points.append(pt)

        return circ_points

    def arcCircleSymbol(_pt, _r, _ang):
        x = _pt.getX()
        y = _pt.getY()
        _ang = int((360-_ang)/10)
        num = 36
        arc_points = []

        for i in range(0, num-_ang+1):
            theta = 2*CompGeom.PI*i/num
            pt = Point(x + _r*math.cos(theta), y + _r*math.sin(theta))
            arc_points.append(pt)

        return arc_points

    @ staticmethod
    def arrowSymbol(_pt, _scale, _ang):

        _ang = _ang*CompGeom.PI/180
        x = AttribSymbols.rotateCoord(Point(3*_scale, 0), _ang)
        y = AttribSymbols.rotateCoord(Point(0, 3*_scale), _ang)

        pt = _pt + x*0.1
        pt_a = pt + x*0.1
        pt_b = pt_a + y*0.1
        pt_c = pt_a - y*0.1
        pt_d = pt + x

        return [pt, pt_d], [pt, pt_b, pt_c]

    def arrowPointCL(_attribute, _pt, _scale):

        properties = _attribute['properties']
        lines = []
        triangles = []
        circles = []

        if properties['Fx'] != 0:
            if properties['Fx'] < 0:
                line, tr = AttribSymbols.arrowSymbol(_pt, _scale, 0)
            else:
                line, tr = AttribSymbols.arrowSymbol(_pt, _scale, 180)
            lines.append(line)
            triangles.append(tr)

        if properties['Fy'] != 0:
            if properties['Fy'] < 0:
                line, tr = AttribSymbols.arrowSymbol(_pt, _scale, 270)

            else:
                line, tr = AttribSymbols.arrowSymbol(_pt, _scale, 90)

            lines.append(line)
            triangles.append(tr)

        if properties['Mz'] != 0:
            if properties['Mz'] < 0:
                cr = AttribSymbols.arcCircleSymbol(_pt, _scale, 180)
                tr = AttribSymbols.triangleSymbol(
                    cr[0], _scale*0.5, 90)
                circles.append(cr)
                triangles.append(tr)
            else:
                cr = AttribSymbols.arcCircleSymbol(_pt, _scale, 180)
                tr = AttribSymbols.triangleSymbol(
                    cr[-1], _scale*0.5, 90)
                circles.append(cr)
                triangles.append(tr)

        return lines, triangles, circles

    def arrowSegmentUL(_attribute, _seg, _scale):
        properties = _attribute['properties']
        lines = []
        triangles = []
        disp = Point(0, 0)
        points = _seg.getPoints().copy()

        while len(points) >= 2:
            aux_line = Line(points[0], points[1])

            if properties['Direction']["index"] == 1:
                local = True
                v = points[1] - points[0]
            else:
                local = False

            if properties['Qx'] != 0:
                if properties['Qx'] > 0:
                    ang = 180
                    if local:
                        ang = AttribSymbols.getAngWithXDirec(v)
                        if points[1].getY() < points[0].getY():
                            ang = ang + 180
                        else:
                            ang = - ang + 180

                    l, tr = AttribSymbols.arrowSegment(
                        _scale*0.45, disp, 0.2, 0.1, 0.9, aux_line, ang, True)
                else:
                    ang = 0
                    if local:
                        ang = AttribSymbols.getAngWithXDirec(v)
                        if not points[1].getY() < points[0].getY():
                            ang = -ang

                    l, tr = AttribSymbols.arrowSegment(
                        _scale*0.45, disp, 0.2, 0.1, 0.9, aux_line, ang, False)

                lines.extend(l)
                triangles.extend(tr)

            if properties['Qy'] != 0:
                if properties['Qy'] > 0:
                    ang = 90
                    if local:
                        ang = AttribSymbols.getAngWithXDirec(v)
                        if points[1].getY() < points[0].getY():
                            ang = ang + 90
                        else:
                            ang = -ang + 90

                    l, tr = AttribSymbols.arrowSegment(
                        _scale*0.5, disp, 0.2, 0, 1, aux_line, ang, True)
                else:
                    ang = 270
                    if local:
                        ang = AttribSymbols.getAngWithXDirec(v)
                        if points[1].getY() < points[0].getY():
                            ang = ang + 270
                        else:
                            ang = -ang + 270

                    l, tr = AttribSymbols.arrowSegment(
                        _scale*0.5, disp, 0.2, 0, 1, aux_line, ang, False)

                lines.extend(l)
                triangles.extend(tr)

            points.pop(0)

        return lines, triangles

    def arrowSegment(_scale, _displc, _step, _init, _end, _seg, _ang, _orient):

        lines = []
        triangles = []
        step = _step
        cont = _init

        if _orient:
            displc = _displc
        else:
            displc = _displc*(-1)

        while cont <= _end:
            pt = _seg.getPoint(cont)
            pt = pt - displc*0.2
            l, tr = AttribSymbols.arrowSymbol(
                pt, _scale, _ang)
            lines.append(l)
            triangles.append(tr)
            cont = cont + step

        return lines, triangles

    @ staticmethod
    def supportPoint(_attribute, _pt, _scale):
        _scale = _scale*0.6
        properties = _attribute['properties']
        x = Point(1*_scale, 0)
        y = Point(0, 1*_scale)
        lines = []
        triangles = []
        squares = []
        circles = []

        if properties['Dx']:

            pt = Point(_pt.getX(), _pt.getY())
            displac = Point(x.getX(), x.getY())

            if properties['Dx pos']["index"] == 0:
                displac = displac*(-1)
                # Left
                if properties['Rz']:
                    pt = pt + displac/4

                tr = AttribSymbols.triangleSymbol(pt, _scale, 0)
                pt_d = tr[1] - x*0.1
                pt_e = tr[2] - x*0.1
            else:
                # Right
                if properties['Rz']:
                    pt = pt + displac/4

                tr = AttribSymbols.triangleSymbol(pt, _scale, 180)
                pt_d = tr[1] + x*0.1
                pt_e = tr[2] + x*0.1

            lines.append([pt_d, pt_e])
            triangles.append(tr)

            if properties['Dx value'] != 0:

                if properties['Dx value'] < 0:

                    if displac.getX() < 0:
                        pt_arrow = (pt_d+pt_e)/2 + displac*2
                    else:
                        pt_arrow = (pt_d+pt_e)/2 + displac/4

                    l, tr = AttribSymbols.arrowSymbol(
                        pt_arrow, _scale*0.5, 0)
                else:
                    if displac.getX() < 0:
                        pt_arrow = (pt_d+pt_e)/2 + displac/4
                    else:
                        pt_arrow = (pt_d+pt_e)/2 + displac*2

                    l, tr = AttribSymbols.arrowSymbol(
                        pt_arrow, _scale*0.5, 180)

                lines.append(l)
                triangles.append(tr)

        if properties['Dy']:

            pt = Point(_pt.getX(), _pt.getY())
            displac = Point(y.getX(), y.getY())

            if properties['Dy pos']["index"] == 0:
                displac = displac*(-1)
                # Down
                if properties['Rz']:
                    pt = pt + displac/4

                tr = AttribSymbols.triangleSymbol(pt, _scale, 270)
                pt_d = tr[1] - y*0.1
                pt_e = tr[2] - y*0.1
            else:
                # Up
                if properties['Rz']:
                    pt = pt + displac/4

                tr = AttribSymbols.triangleSymbol(pt, _scale, 90)
                pt_d = tr[1] + y*0.1
                pt_e = tr[2] + y*0.1

            lines.append([pt_d, pt_e])
            triangles.append(tr)

            if properties['Dy value'] != 0:

                if properties['Dy value'] < 0:

                    if displac.getY() < 0:
                        pt_arrow = (pt_d+pt_e)/2 + displac*2
                    else:
                        pt_arrow = (pt_d+pt_e)/2 + displac/4

                    l, tr = AttribSymbols.arrowSymbol(
                        pt_arrow, _scale*0.5, 270)
                else:
                    if displac.getY() < 0:
                        pt_arrow = (pt_d+pt_e)/2 + displac/4
                    else:
                        pt_arrow = (pt_d+pt_e)/2 + displac*2

                    l, tr = AttribSymbols.arrowSymbol(
                        pt_arrow, _scale*0.5, 90)

                lines.append(l)
                triangles.append(tr)

        if properties['Rz']:
            sq = AttribSymbols.squareSymbol(_pt, _scale, 0)
            squares.append(sq)

            if properties['Rz value'] != 0:
                if properties['Rz value'] < 0:
                    cr = AttribSymbols.arcCircleSymbol(_pt, _scale*1.4, 180)
                    tr = AttribSymbols.triangleSymbol(
                        cr[0], _scale*0.5, 90)

                else:
                    cr = AttribSymbols.arcCircleSymbol(_pt, _scale*1.4, 180)
                    tr = AttribSymbols.triangleSymbol(
                        cr[-1], _scale*0.5, 90)

                circles.append(cr)
                triangles.append(tr)

        return lines, triangles, squares, circles

    @ staticmethod
    def supportSegment(_attribute, _seg, _scale):

        lines = []
        triangles = []
        squares = []
        circles = []
        seg_pts = _seg.getPoints()
        points = []
        points.append(seg_pts[0])
        points.append(seg_pts[-1])
        points.append(_seg.getPoint(0.5))

        for pt in points:
            l, tr, sq, circ = AttribSymbols.supportPoint(
                _attribute, pt, _scale)
            lines.extend(l)
            triangles.extend(tr)
            squares.extend(sq)
            circles.extend(circ)

        return lines, triangles, squares, circles

    def Nsbdvs(_attribute, _seg):
        points = []
        properties = _attribute['properties']
        nsudv = properties['Value']
        ratio = properties['Ratio']
        points = CompGeom.getNumberOfSudvisions(_seg, nsudv, ratio, False)

        return points


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ------------------------ COMPUTACIONAL GEOMETRY ----------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


class CompGeom:

    ABSTOL = 1e-7  # Absolute tolerance value
    PI = 3.1415926535897932384626433832975028841971693993751058

    @staticmethod
    def orient2d(pa, pb, pc):
        acx = pa[0] - pc[0]
        bcx = pb[0] - pc[0]
        acy = pa[1] - pc[1]
        bcy = pb[1] - pc[1]
        return acx * bcy - acy * bcx

    # Return the sign (NEGATIVE, ZERO, or POSITIVE) of the oriented
    # twice area formed by three given points.
    @staticmethod
    def signOrient2d(_p1, _p2, _p3):

        det = 0.0
        pa = [_p1.getX(), _p1.getY()]
        pb = [_p2.getX(), _p2.getY()]
        pc = [_p3.getX(), _p3.getY()]

        det = CompGeom.orient2d(pa, pb, pc)
        if(det != 0.0):
            if(det > 0.0):
                return "POSITIVE"
            else:
                return "NEGATIVE"

        return "ZERO"

    # Return the signed value of the oriented twice area formed by
    # three given points.
    @staticmethod
    def valOrient2d(_p1, _p2, _p3):
        pa = [_p1.getX(), _p1.getY()]
        pb = [_p2.getX(), _p2.getY()]
        pc = [_p3.getX(), _p3.getY()]

        det = CompGeom.orient2d(pa, pb, pc)
        return det

    # Return a flag (true or false) stating whether the three given
    # points are collinear.
    @staticmethod
    def areCollinear(_p1,  _p2, _p3):
        return CompGeom.signOrient2d(_p1, _p2, _p3) == 'ZERO'

    # Return a flag (true or false) stating whether the third given
    # point is on the left side of oriented segment formed by the
    # first given points.
    @staticmethod
    def isLeftSide(_p1, _p2, _p3):
        return CompGeom.signOrient2d(_p1, _p2, _p3) == 'POSITIVE'

    # Return a flag (true or false) stating whether the third given
    # point is on the right side of oriented segment formed by the
    # first given points.
    @staticmethod
    def isRightSide(_p1, _p2, _p3):
        return CompGeom.signOrient2d(_p1, _p2, _p3) == 'NEGATIVE'

    # Return the sign (NEGATIVE, ZERO, or POSITIVE) of the oriented
    # twice area formed by three given points.
    # This function uses conventional floating-point operations to
    # compute the area and compares the area result with a hard-coded
    # very small tolerance value (ABSTOL).
    @staticmethod
    def signArea2d(_p1, _p2, _p3):
        det = Point.area2d(_p1, _p2, _p3)
        if abs(det) < CompGeom.ABSTOL:
            return "ZERO"
        if det > 0.0:
            return "POSITIVE"
        return "NEGATIVE"

    # Return the signed value of the oriented twice area formed by
    # three given points.
    # This function uses conventional floating-point operations to
    # compute the signed area.
    @staticmethod
    def valArea2d(_p1,  _p2, _p3):
        return Point.area2d(_p1, _p2, _p3)

    # Get closest point on line 'p1'-'p2'.
    # Returns the distance between given point and closest point.
    # Also returns parametric value (from -infinity to +infinity)
    # of closest point along the line.
    @staticmethod
    def getClosestPointLine(_p1, _p2, _p, _pC, _t):
        v12 = _p2 - _p1
        v1p = _p - _p1

        _t = Point.dotprod(v12, v1p) / Point.sizesquare(v12)
        _pC = _p1 + v12 * _t

        dist = Point.size(_p - _pC)
        return dist

    # Get closest point on line segment 'p1'-'p2'.
    # Returns the distance between given point and closest point.
    # Also returns parametric value (between 0 and 1) of closest point
    # along the line.
    # The difference between this function and function 'getClosestPointLine'
    # is that a closest point outside the limits of segment 'p1-p2' is
    # snapped to one of the segment end points.
    @staticmethod
    def getClosestPointSegment(_p1, _p2, _p):
        v12 = _p2 - _p1
        v1p = _p - _p1

        t = Point.dotprod(v12, v1p) / Point.sizesquare(v12)

        if abs(t) < CompGeom.ABSTOL or t < 0.0:
            pC = _p1
            t = 0.0
        elif abs(t-1.0) < CompGeom.ABSTOL or t > 1.0:
            pC = _p2
            t = 1.0
        else:
            pC = _p1 + v12 * t

        dist = Point.euclidiandistance(pC, _p)
        return dist, pC, t

    # Check for collinear segments 'p1'-'p2' and 'p3'-'p4'.
    @staticmethod
    def checkCollinearSegments(_p1, _p2, _p3, _p4):
        sign123 = CompGeom.signArea2d(_p1, _p2, _p3)
        sign124 = CompGeom.signArea2d(_p1, _p2, _p4)

        # Check for collinear segments
        if sign123 == 'ZERO' and sign124 == 'ZERO':
            return True
        return False

    # Checks for two line segments intersection:
    # Checks whether segment 'p1'-'p2' intercepts segment 'p3'-'p4'.
    # Returns an integer intersection type value.
    # In case there is intersection, outputs the result in 'pi' parameter
    # and returns parametric values ('t12' and 't34' between 0 and 1) along
    # the two segments.
    # Ref.:
    # M. Gavrilova & J.G. Rokne - Reliable line segment intersection testing,
    # Computer-Aided Design, Vol. 32, Issue 12, pp. 737�745, 2000.
    @staticmethod
    def computeSegmentSegmentIntersection(_p1, _p2, _p3, _p4):
        # Discard intersection if second segment is located to the left (_l) or
        # to the right (_r) of horizontal bounding box limits of first segment.
        x12_l = min(_p1.getX(), _p2.getX())
        x12_r = max(_p1.getX(), _p2.getX())
        x34_l = min(_p3.getX(), _p4.getX())
        x34_r = max(_p3.getX(), _p4.getX())

        if x12_r+CompGeom.ABSTOL < x34_l or x34_r < x12_l-CompGeom.ABSTOL:
            return "DO_NOT_INTERSECT", None, None, None

        # Discard intersection if second segment is located below
        #  bottom (_b) or above top (_t) of bounding box of first segment.
        y12_b = min(_p1.getY(), _p2.getY())
        y12_t = max(_p1.getY(), _p2.getY())
        y34_b = min(_p3.getY(), _p4.getY())
        y34_t = max(_p3.getY(), _p4.getY())

        if y12_t + CompGeom.ABSTOL < y34_b or y34_t < y12_b-CompGeom.ABSTOL:
            return 'DO_NOT_INTERSECT', None, None, None

        # Get signs of oriented twice area for points p1-p2-p3 and for points p1-p2-p4
        sign123 = CompGeom.signArea2d(_p1, _p2, _p3)
        sign124 = CompGeom.signArea2d(_p1, _p2, _p4)

        # Check for collinear segments
        if sign123 == 'ZERO' and sign124 == 'ZERO':
            return 'COLLINEAR', None, None, None

        # Check for second segment on the same side of first segment
        if ((sign123 == 'POSITIVE' and sign124 == 'POSITIVE') or
                (sign123 == 'NEGATIVE' and sign124 == 'NEGATIVE')):
            return 'DO_NOT_INTERSECT', None, None, None

        # Get signs of oriented twice area for points p3-p4-p1 and for points p3-p4-p2
        sign341 = CompGeom.signArea2d(_p3, _p4, _p1)
        sign342 = CompGeom.signArea2d(_p3, _p4, _p2)

        # Check for first segment on the same side of second segment
        if ((sign341 == 'POSITIVE') and (sign342 == 'POSITIVE') or
                (sign341 == 'NEGATIVE' and sign342 == 'NEGATIVE')):
            return 'DO_NOT_INTERSECT', None, None, None

        # Check for one point of second segment touching first segment.
        # Also compute the intersection point and the parametric values
        # ('t12' and 't34' between 0 and 1) along the two segments.
        # In this case, 't34' is either equal to 0 or equal to 1.
        if sign123 == 'ZERO' or sign124 == 'ZERO':
            if sign123 == 'ZERO':
                t34 = 0.0
                pi = _p3
            elif sign124 == 'ZERO':
                t34 = 1.0
                pi = _p4

            if sign341 == 'ZERO':
                t12 = 0.0
                pi = _p1
            elif sign342 == 'ZERO':
                t12 = 1.0
                pi = _p2
            else:
                area341 = CompGeom.valArea2d(_p3, _p4, _p1)
                area342 = CompGeom.valArea2d(_p3, _p4, _p2)
                t12 = area341 / (area341 - area342)

            return 'TOUCH', pi, t12, t34

        # Check for one point of first segment touching second segment
        # Also compute the intersection point and the parametric values
        # ('t12' and 't34' between 0 and 1) along the two segments.
        # In this case, 't12' is either equal to 0 or equal to 1.
        if sign341 == 'ZERO' or sign342 == 'ZERO':
            if sign341 == 'ZERO':
                t12 = 0.0
                pi = _p1
            elif sign342 == 'ZERO':
                t12 = 1.0
                pi = _p2

            if sign123 == 'ZERO':
                t34 = 0.0
                pi = _p3
            elif sign124 == 'ZERO':
                t34 = 1.0
                pi = _p4
            else:
                area123 = CompGeom.valArea2d(_p1, _p2, _p3)
                area124 = CompGeom.valArea2d(_p1, _p2, _p4)
                t34 = area123 / (area123 - area124)

            return 'TOUCH', pi, t12, t34

        # When get to this point, there is an intersection point of the
        # two segments. Compute parametric values of intersection point
        # along each segment.
        area341 = CompGeom.valArea2d(_p3, _p4, _p1)
        area342 = CompGeom.valArea2d(_p3, _p4, _p2)
        area123 = CompGeom.valArea2d(_p1, _p2, _p3)
        area124 = CompGeom.valArea2d(_p1, _p2, _p4)
        t12 = area341 / (area341 - area342)
        t34 = area123 / (area123 - area124)

        # Compute intersection point (there are two equivalent options)
        v34 = _p4 - _p3
        _pi = _p3 + v34*t34
        return 'DO_INTERSECT', _pi, t12, t34

    # This function classifies the projection of a given 'p' point w.r.t. a
    # given segment 'p1-p2'.
    # It returns the classified position of the projection of a point
    # along the infinite line that contains a given segment.
    # It also returns the parametric value of the project point along segment.
    # Its main used is to classify the points of two collinear segments,
    # but it may be used to classify points projected at the infinite line
    # that contains a segment.
    @staticmethod
    def getPtPosWrtSegment(_p1, _p2, _p):
        v12 = _p2 - _p1
        v1p = _p - _p1
        # Get parametric value of project point on segment line
        _t = Point.dotprod(v12, v1p) / Point.sizesquare(v12)

        if abs(_t) < CompGeom.ABSTOL:
            return 'START_SEG', _t  # At start point of segment
        elif abs(_t-1.0) < CompGeom.ABSTOL:
            return 'END_SEG', _t   # At end point of segment
        elif _t < 0.0:
            return 'BEFORE_SEG', _t   # Outside and before segment
        elif _t > 1.0:
            return 'AFTER_SEG', _t   # Outside and after segment
        return 'INSIDE_SEG', _t   # Inside segment

    # This function returns a flag indicating whether the vertices
    # of a polygon are in counter-clockwise order.
    # The algorithm is as follows:
    # Traverse the loop of coordinates, assuming that it is in counter-
    # clockwise order, computing the components of the "area" of the
    # enclosed polygon.  The total "area" components are computed by
    # adding "area" components (cross product components) of
    # triangles sides formed by the first, previous, and current
    # vertices.  If the loop is not convex, some of the triangle
    # areas might be negative, but those will be compensated by other
    # positive triangle areas so that the final area is positive.
    # (area here is actually twice the area).
    # (positive here means in the direction of the face normal).
    @staticmethod
    def isCounterClockwisePolygon(_poly):
        area = 0.0  # twice the enclosed polygon area

        # Compute area assuming that polygon is in counter-clockwise order
        for i in range(2, len(_poly)):
            area += Point.area2d(_poly[0], _poly[i-1], _poly[i])

        # If area is greater than zero, then polygon is in counter-clockwise
        # order, otherwise it is in clockwise order.
        if area > 0.0:
            return True
        return False

    # This function returns a flag indicating whether the given point 'p'
    # is inside a polygon 'poly'.
    # The algorithm counts the number of intersections that a horizontal
    # line emanating at the point 'p' in positive x direction makes with
    # the boundary lines of the polygon. If the number of intersections
    # is odd, the point is inside de polygon. Otherwise, it is outside
    # the polygon.
    # If ray passes at an horizontal segment on polygon boundary, do not
    # count any intersection.
    # If ray passes at a vertex of the polygon, only consider intersection
    # if boundary segment is above ray.
    @staticmethod
    def isPointInPolygon(_poly, _p):
        x = _p.getX()
        y = _p.getY()
        n = len(_poly)  # number of polygon points
        ni = 0  # number of intersections

        for i in range(0, n):
            p1 = _poly[i]  # first point of current line segment
            p2 = _poly[(i+1) % n]  # second point of current line segment

            if (p1.getY() == p2.getY()):  # discard horizontal line
                continue

            if p1.getY() > y and p2.getY() > y:  # discard line above ray
                continue

            if p1.getX() < x and p2.getX() < x:  # discard line to the left of point
                continue

            if p1.getY() < y and p2.getY() < y:  # Discard line below ray
                continue

            if p1.getY() == y:  # ray passes at first line point
                if p1.getX() > x and p2.getY() > y:
                    # Count intersection if first point is to the right of given point
                    # and second point is above.
                    ni += 1
            else:
                if p2.getY() == y:  # ray passes at second point
                    if p2.getX() > x and p1.getY() > y:
                        # Count intersection if first point is to the right of given point
                        # and second point is above.
                        ni += 1
                else:  # ray passes with first and second points
                    if p1.getX() > x and p2.getX() > x:
                        # Count intersection if first point is to the right of given point
                        # and second point is above.
                        ni += 1
                    else:
                        # Compute x coordinate of intersection of ray with line segment
                        dx = p1.getX() - p2.getX()
                        xc = p1.getX()

                        if dx != 0.0:
                            xc += (y - p1.getY())*dx / (p1.getY()-p2.getY())

                        if xc > x:
                            # Count intersection if first point is to the right of given point
                            # and second point is above.
                            ni += 1

        # If number of intersections is odd, point is inside polygon.
        if (ni % 2) > 0:
            return True

        # If number of intersections if even, point is outside polygon.
        return False

    @staticmethod
    def computeLineIntersection(_p1, _p2, _p3, _p4):

        pts = []
        param1 = []
        param2 = []

        status, pi, t12, t34 = CompGeom.computeSegmentSegmentIntersection(
            _p1, _p2, _p3, _p4)

        if status == 'DO_NOT_INTERSECT':
            return False, pts, param1, param2

        elif status == 'DO_INTERSECT':
            pts.append(pi)
            param1.append(t12)
            param2.append(t34)
            return True, pts, param1, param2

        elif status == 'COLLINEAR':
            pos3_12, t3_12 = CompGeom.getPtPosWrtSegment(_p1, _p2, _p3)
            pos4_12, t4_12 = CompGeom.getPtPosWrtSegment(_p1, _p2, _p4)
            pos1_34, t1_34 = CompGeom.getPtPosWrtSegment(_p3, _p4, _p1)
            pos2_34, t2_34 = CompGeom.getPtPosWrtSegment(_p3, _p4, _p2)

            if ((pos3_12 == 'BEFORE_SEG' and pos4_12 == 'BEFORE_SEG') or
                    (pos3_12 == 'AFTER_SEG' and pos4_12 == 'AFTER_SEG')):
                # The two segments do not intercept
                return False, pts, param1, param2

            elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'START_SEG':
                pts.append(_p1)
                param1.append(0.0)
                param2.append(1.0)
                return True, pts, param1, param2

            elif pos3_12 == 'START_SEG' and pos4_12 == 'BEFORE_SEG':
                pts.append(_p1)
                param1.append(0.0)
                param2.append(0.0)
                return True, pts, param1, param2

            elif pos3_12 == 'END_SEG' and pos4_12 == 'AFTER_SEG':
                pts.append(_p2)
                param1.append(1.0)
                param2.append(0.0)
                return True, pts, param1, param2

            elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'END_SEG':
                pts.append(_p2)
                param1.append(1.0)
                param2.append(1.0)
                return True, pts, param1, param2

            elif pos3_12 == 'START_SEG' and pos4_12 == 'END_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(0.0)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(1.0)
                return True, pts, param1, param2

            elif pos3_12 == 'END_SEG1' and pos4_12 == 'START_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(1.0)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(0.0)
                return True, pts, param1, param2

            elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'INSIDE_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(t1_34)

                # Store second pair of intersection parameters
                pts.append(_p4)
                param1.append(t4_12)
                param2.append(1.0)
                return True, pts, param1, param2

            elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'BEFORE_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(t1_34)

                # Store second pair of intersection parameters
                pts.append(_p3)
                param1.append(t3_12)
                param2.append(0.0)
                return True, pts, param1, param2

            elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'END_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(t1_34)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(1.0)
                return True, pts, param1, param2

            elif pos3_12 == 'END_SEG' and pos4_12 == 'BEFORE_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(t1_34)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(0.0)
                return True, pts, param1, param2

            elif pos3_12 == 'START_SEG' and pos4_12 == 'INSIDE_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(0.0)

                # Store second pair of intersection parameters
                pts.append(_p4)
                param1.append(t4_12)
                param2.append(1.0)
                return True, pts, param1, param2

            elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'START_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(1.0)

                # Store second pair of intersection parameters
                pts.append(_p3)
                param1.append(t3_12)
                param2.append(0.0)
                return True, pts, param1, param2

            elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'INSIDE_SEG':

                if t3_12 < t4_12:
                    # Store fisrt pair of intersection parameters
                    pts.append(_p3)
                    param1.append(t3_12)
                    param2.append(0.0)

                    # Store second pair of intersection parameters
                    pts.append(_p4)
                    param1.append(t4_12)
                    param2.append(1.0)

                else:
                    # Store fisrt pair of intersection parameters
                    pts.append(_p4)
                    param1.append(t4_12)
                    param2.append(1.0)

                    # Store second pair of intersection parameters
                    pts.append(_p3)
                    param1.append(t3_12)
                    param2.append(0.0)

                return True, pts, param1, param2

            elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'AFTER_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(t1_34)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(t2_34)
                return True, pts, param1, param2

            elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'BEFORE_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(t1_34)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(t2_34)
                return True, pts, param1, param2

            elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'END_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p3)
                param1.append(t3_12)
                param2.append(0.0)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(1.0)
                return True, pts, param1, param2

            elif pos3_12 == 'END_SEG' and pos4_12 == 'INSIDE_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p4)
                param1.append(t4_12)
                param2.append(1.0)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(0.0)
                return True, pts, param1, param2

            elif pos3_12 == 'START_SEG' and pos4_12 == 'AFTER_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(0.0)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(t2_34)
                return True, pts, param1, param2

            elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'START_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p1)
                param1.append(0.0)
                param2.append(1.0)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(t2_34)
                return True, pts, param1, param2

            elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'AFTER_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p3)
                param1.append(t3_12)
                param2.append(0.0)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(t2_34)
                return True, pts, param1, param2

            elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'INSIDE_SEG':
                # Store fisrt pair of intersection parameters
                pts.append(_p4)
                param1.append(t4_12)
                param2.append(1.0)

                # Store second pair of intersection parameters
                pts.append(_p2)
                param1.append(1.0)
                param2.append(t2_34)
                return True, pts, param1, param2

        elif status == 'TOUCH':
            # one segments touches the other, in the middle or extremity!
            # checks if the curves are not consecutive segments
            # checks if the polygons touch at the extremities
            pts.append(pi)
            param1.append(t12)
            param2.append(t34)
            return True, pts, param1, param2

        return False, pts, param1, param2

    @staticmethod
    def splitSelfIntersected(_poly):

        # verifies for each pair of possible segments if they intersect, and
        # stores for both segments the parametric coordinate where intersection occurs
        iStatus = False
        segONETotalLength = 0.0
        segTWOTotalLength = 0.0
        intersecParams = []
        params = []
        pts = []
        for i in range(0, len(_poly)-1):
            segONEPartialLength = Point.euclidiandistance(
                _poly[i], _poly[i + 1])
            segTWOTotalLength = segONETotalLength + segONEPartialLength
            segTWOPartialLength = 0.0

            for j in range(i+1, len(_poly)-1):
                segTWOPartialLength = Point.euclidiandistance(
                    _poly[j], _poly[j + 1])
                status, pi, t12, t34 = CompGeom.computeSegmentSegmentIntersection(_poly[i], _poly[i + 1],
                                                                                  _poly[j], _poly[j + 1])

                if status == 'DO_NOT_INTERSECT':
                    # do nothing, continue the checking!
                    pass
                elif status == 'DO_INTERSECT':
                    # the straight segments intersect in the middle!
                    intersecParams.append([
                        segONETotalLength + t12*segONEPartialLength, pi])
                    intersecParams.append([
                        segTWOTotalLength + t34*segTWOPartialLength, pi])
                    iStatus = True
                elif status == 'COLLINEAR':
                    # the straight segments are collinear !
                    pos3_12, t3_12 = CompGeom.getPtPosWrtSegment(
                        _poly[i], _poly[i + 1], _poly[j])
                    pos4_12, t4_12 = CompGeom.getPtPosWrtSegment(
                        _poly[i], _poly[i + 1], _poly[j + 1])
                    pos1_34, t1_34 = CompGeom.getPtPosWrtSegment(
                        _poly[j], _poly[j + 1], _poly[i])
                    pos2_34, t2_34 = CompGeom.getPtPosWrtSegment(
                        _poly[j], _poly[j + 1], _poly[i + 1])

                    if(pos3_12 == 'BEFORE_SEG' and pos4_12 == 'BEFORE_SEG' or
                       pos3_12 == 'AFTER_SEG' and pos4_12 == 'AFTER_SEG'):
                        # The two segments do not intercept
                        pass
                    elif((pos3_12 == 'BEFORE_SEG' and pos4_12 == 'START_SEG') or
                         (pos3_12 == 'START_SEG' and pos4_12 == 'BEFORE_SEG') or
                         (pos3_12 == 'END_SEG' and pos4_12 == 'AFTER_SEG') or
                         (pos3_12 == 'AFTER_SEG' and pos4_12 == 'END_SEG')):

                        # Segments simply touch at one end without overlapping
                        if i == 0 and j == len(_poly)-2:
                            segONEInterAtParam = segONETotalLength
                            segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength
                            intersecParams.append(
                                [segONEInterAtParam, _poly[i]])
                            intersecParams.append(
                                [segTWOInterAtParam, _poly[j + 1]])
                            iStatus = True

                    elif pos3_12 == 'START_SEG' and pos4_12 == 'END_SEG':
                        # Segments have common end points: just delete second segment.
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, _poly[i + 1]])
                        intersecParams.append(
                            [segTWOInterAtParam, _poly[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'END_SEG' and pos4_12 == 'START_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength
                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j]])
                        intersecParams.append([segTWOInterAtParam, _poly[j]])
                        iStatus = True

                    elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'INSIDE_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[j+1]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'BEFORE_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j]])
                        intersecParams.append([segTWOInterAtParam, _poly[j]])
                        iStatus = True

                    elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'END_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[i+1]])
                        iStatus = True

                    elif pos3_12 == 'END_SEG' and pos4_12 == 'BEFORE_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j]])
                        intersecParams.append([segTWOInterAtParam, _poly[j]])
                        iStatus = True

                    elif pos3_12 == 'START_SEG' and pos4_12 == 'INSIDE_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[j+1]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'START_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j]])
                        intersecParams.append([segTWOInterAtParam, _poly[j]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'INSIDE_SEG':
                        if t3_12 < t4_12:
                            segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                            segTWOInterAtParam = segTWOTotalLength

                            # Store fisrt pair of intersection parameters
                            intersecParams.append(
                                [segONEInterAtParam, _poly[j]])
                            intersecParams.append(
                                [segTWOInterAtParam, _poly[j]])

                            segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                            segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                            # Store second pair of intersection parameters
                            intersecParams.append(
                                [segONEInterAtParam, _poly[j+1]])
                            intersecParams.append(
                                [segTWOInterAtParam, _poly[j+1]])

                        else:
                            segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                            segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                            # Store fisrt pair of intersection parameters
                            intersecParams.append(
                                [segONEInterAtParam, _poly[j+1]])
                            intersecParams.append(
                                [segTWOInterAtParam, _poly[j+1]])

                            segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                            segTWOInterAtParam = segTWOTotalLength

                            # Store second pair of intersection parameters
                            intersecParams.append(
                                [segONEInterAtParam, _poly[j]])
                            intersecParams.append(
                                [segTWOInterAtParam, _poly[j]])

                        iStatus = True

                    elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'AFTER_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[i+1]])
                        iStatus = True

                    elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'BEFORE_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[i+1]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'END_SEG':
                        segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j]])
                        intersecParams.append([segTWOInterAtParam, _poly[j]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[i+1]])
                        iStatus = True

                    elif pos3_12 == 'END_SEG' and pos4_12 == 'INSIDE_SEG':
                        segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[j+1]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j]])
                        intersecParams.append([segTWOInterAtParam, _poly[j]])
                        iStatus = True

                    elif pos3_12 == 'START_SEG' and pos4_12 == 'AFTER_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[i+1]])
                        iStatus = True

                    elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'START_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i]])
                        intersecParams.append([segTWOInterAtParam, _poly[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[i+1]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'AFTER_SEG':
                        segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j]])
                        intersecParams.append([segTWOInterAtParam, _poly[j]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[i+1]])
                        iStatus = True

                    elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'INSIDE_SEG':
                        segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[j+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[j+1]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append([segONEInterAtParam, _poly[i+1]])
                        intersecParams.append([segTWOInterAtParam, _poly[i+1]])
                        iStatus = True
                elif status == 'TOUCH':
                    # one segments touches the other

                    # avoid consecutive segments
                    if j != (i+1):
                        intersecParams.append(
                            [segONETotalLength + t12*segONEPartialLength, pi])
                        intersecParams.append(
                            [segTWOTotalLength + t34*segTWOPartialLength, pi])
                        iStatus = True

                segTWOTotalLength += segTWOPartialLength

            segONETotalLength += segONEPartialLength

        # removes duplicate elements
        unique_intersecParams = []
        for item in intersecParams:
            if item not in unique_intersecParams:
                unique_intersecParams.append(item)

        unique_intersecParams.sort()

        # Calculate the parameters based on its partial length
        for it in unique_intersecParams:
            params.append(it[0]/segONETotalLength)
            pts.append(it[1])

        return iStatus, pts, params

    @staticmethod
    def computePolyPolyIntersection(_poly1, _poly2):

        # verifies for each pair of possible segments if they intersect, and
        # stores for both segments the parametric coordinate where intersection occurs
        segONETotalLength = 0.0
        iStatus = False
        intersecParams = []
        param1 = []
        param2 = []
        pts = []

        for i in range(0, len(_poly1)-1):
            segONEPartialLength = Point.euclidiandistance(
                _poly1[i], _poly1[i + 1])
            segTWOPartialLength = 0.0
            segTWOTotalLength = 0.0

            for j in range(0, len(_poly2)-1):

                segTWOPartialLength = Point.euclidiandistance(
                    _poly2[j], _poly2[j + 1])
                status, pi, t12, t34 = CompGeom.computeSegmentSegmentIntersection(
                    _poly1[i], _poly1[i+1], _poly2[j], _poly2[j + 1])

                if status == 'DO_NOT_INTERSECT':
                    # do nothing, continue the checking
                    pass
                elif status == 'DO_INTERSECT':
                    # the straight segments intersect in the middle!
                    segONEInterAtParam = segONETotalLength + t12*segONEPartialLength
                    segTWOInterAtParam = segTWOTotalLength + t34*segTWOPartialLength
                    intersecParams.append(
                        [segONEInterAtParam, segTWOInterAtParam, pi])
                    iStatus = True

                elif status == 'COLLINEAR':
                    pos3_12, t3_12 = CompGeom.getPtPosWrtSegment(
                        _poly1[i], _poly1[i + 1], _poly2[j])
                    pos4_12, t4_12 = CompGeom.getPtPosWrtSegment(
                        _poly1[i], _poly1[i + 1], _poly2[j + 1])
                    pos1_34, t1_34 = CompGeom.getPtPosWrtSegment(
                        _poly2[j], _poly2[j + 1], _poly1[i])
                    pos2_34, t2_34 = CompGeom.getPtPosWrtSegment(
                        _poly2[j], _poly2[j + 1], _poly1[i + 1])

                    if ((pos3_12 == 'BEFORE_SEG' and pos4_12 == 'BEFORE_SEG') or
                            (pos3_12 == 'AFTER_SEG' and pos4_12 == 'AFTER_SEG')):
                        # The two segments do not intercept
                        pass

                    elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'START_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])
                        iStatus = True

                    elif pos3_12 == 'START_SEG' and pos4_12 == 'BEFORE_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])
                        iStatus = True

                    elif pos3_12 == 'END_SEG' and pos4_12 == 'AFTER_SEG':
                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i+1]])
                        iStatus = True

                    elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'END_SEG':
                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'START_SEG' and pos4_12 == 'END_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'END_SEG' and pos4_12 == 'START_SEG':
                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'INSIDE_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly2[j + 1]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'BEFORE_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly2[j]])
                        iStatus = True

                    elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'END_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'END_SEG' and pos4_12 == 'BEFORE_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'START_SEG' and pos4_12 == 'INSIDE_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly2[j + 1]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'START_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly2[j]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'INSIDE_SEG':

                        if t3_12 < t4_12:

                            segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                            segTWOInterAtParam = segTWOTotalLength

                            # Store fisrt pair of intersection parameters
                            intersecParams.append(
                                [segONEInterAtParam, segTWOInterAtParam, _poly2[j]])

                            segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                            segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                            # Store second pair of intersection parameters
                            intersecParams.append(
                                [segONEInterAtParam, segTWOInterAtParam, _poly2[j + 1]])

                        else:

                            segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                            segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                            # Store fisrt pair of intersection parameters
                            intersecParams.append(
                                [segONEInterAtParam, segTWOInterAtParam, _poly2[j + 1]])

                            segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                            segTWOInterAtParam = segTWOTotalLength

                            # Store second pair of intersection parameters
                            intersecParams.append(
                                [segONEInterAtParam, segTWOInterAtParam, _poly2[j]])

                        iStatus = True

                    elif pos3_12 == 'BEFORE_SEG' and pos4_12 == 'AFTER_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'BEFORE_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + t1_34*segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'END_SEG':

                        segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly2[j]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'END_SEG' and pos4_12 == 'INSIDE_SEG':

                        segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly2[j + 1]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'START_SEG' and pos4_12 == 'AFTER_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'START_SEG':

                        segONEInterAtParam = segONETotalLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'INSIDE_SEG' and pos4_12 == 'AFTER_SEG':

                        segONEInterAtParam = segONETotalLength + t3_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly2[j]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                    elif pos3_12 == 'AFTER_SEG' and pos4_12 == 'INSIDE_SEG':

                        segONEInterAtParam = segONETotalLength + t4_12*segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + segTWOPartialLength

                        # Store fisrt pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly2[j + 1]])

                        segONEInterAtParam = segONETotalLength + segONEPartialLength
                        segTWOInterAtParam = segTWOTotalLength + t2_34*segTWOPartialLength

                        # Store second pair of intersection parameters
                        intersecParams.append(
                            [segONEInterAtParam, segTWOInterAtParam, _poly1[i + 1]])
                        iStatus = True

                elif status == 'TOUCH':

                    # one segments touches the other, in the middle or extremity!
                    # checks if the curves are not consecutive segments
                    # checks if the polygons touch at the extremities

                    segONEInterAtParam = segONETotalLength + t12*segONEPartialLength
                    segTWOInterAtParam = segTWOTotalLength + t34*segTWOPartialLength
                    intersecParams.append(
                        [segONEInterAtParam, segTWOInterAtParam, pi])
                    iStatus = True

                segTWOTotalLength += segTWOPartialLength

            segONETotalLength += segONEPartialLength

        # removes duplicate elements
        unique_intersecParams = []
        for item in intersecParams:
            if item not in unique_intersecParams:
                unique_intersecParams.append(item)

        # sorts the pairs of params by the _poly1 parametric order
        unique_intersecParams.sort()

        for it in unique_intersecParams:
            param1.append(it[0]/segONETotalLength)
            param2.append(it[1]/segTWOTotalLength)
            pts.append(it[2])

        return iStatus, pts, param1, param2

    @staticmethod
    def SdvSubdivideSegment(_p1, _p2, _nsudv, _quad, _ratio):
        coords = []
        n_pts = _nsudv-1
        if _quad:
            n_pts = (_nsudv*2) - 1

        for i in range(0, n_pts):
            coords.insert(i, Point())

        # find the edge endpoint coordinates and get the vertices coordinates
        x0 = _p1.x
        y0 = _p1.y
        x1 = _p2.x
        y1 = _p2.y

        _ratio = 1.0 / _ratio
        a = (2.0 * _ratio) / ((_ratio + 1.0) * _nsudv)
        b = (a * (1.0 - _ratio)) / (2.0 * _ratio * (_nsudv - 1.0))

        if _quad:
            # --------------- subdivision for quadratic type elements -------------
            # get the coordinates at the segment ends
            j = 1
            for i in range(1, _nsudv):
                v = a * i + b * i * (i - 1.0)
                u = 1.0 - v
                coords[j].x = u * x0 + v * x1
                coords[j].y = u * y0 + v * y1
                j += 2

            #  --------- get the coordinates at the mid segment points ------------
            # first point on the first segment
            coords[0].x = (x0 + coords[1].x) * 0.5
            coords[0].y = (y0 + coords[1].y) * 0.5

            # now the points on the center segments
            j = 2
            for i in range(1, _nsudv-1):
                coords[j].x = (coords[j-1].x + coords[j+1].x) * 0.5
                coords[j].y = (coords[j-1].y + coords[j+1].y) * 0.5
                j += 2

            # now the point on the last segment
            coords[2*_nsudv-2].x = (coords[2*_nsudv-3].x + x1) * 0.5
            coords[2*_nsudv-2].y = (coords[2*_nsudv-3].y + y1) * 0.5
        else:
            #  ---------- subdivision for non-quadratic type elements ------------
            # get the coordinates at the segment ends
            for i in range(1, _nsudv):
                v = a * i + b * i * (i - 1.0)
                u = 1.0 - v
                coords[i-1].x = u * x0 + v * x1
                coords[i-1].y = u * y0 + v * y1

        return coords

    @staticmethod
    def getNumberOfSudvisions(_segment, _nsbdv, _ratio, _quad):
        length = _segment.length(0, 1)
        coords = []

        if _nsbdv == 0 or _ratio == 0:
            return coords
        elif _nsbdv == 1:
            if _quad:
                coords.append(_segment.getPoint(0.5))
            return coords

        # calculates as if it were straight to find interpolation parameters
        r1 = Point(0.0, 0.0)
        r2 = Point(length, 0.0)
        coords = CompGeom.SdvSubdivideSegment(r1, r2, _nsbdv, _quad, _ratio)

        # calculates the coordinates of the real points on the segment
        pts = []
        for coord in coords:
            t = coord.x / length
            pts.append(_segment.getPoint(t))

        return pts


class Tesselation:
    @staticmethod
    def triangleParing(_p):
        triangs = []
        pn = len(_p)
        left = []  # left neighbor indices
        right = []  # right neighbor indices
        isPolyPt = []

        for i in range(0, pn):
            # initialization
            left.append(((i-1) + pn) % pn)
            right.append(((i+1) + pn) % pn)
            isPolyPt.append(True)

        i = pn-1  # counter

        while len(triangs) < (pn-2):

            i = right[i]
            if (Tesselation.ear_Q(left[i], i, right[i], _p, isPolyPt)):
                # Original implementation (SKIENA & REVILLA, 2002):
                # add_triangle(t,l[i],i,r[i],p);
                tri = [None, None, None]

                tri[0] = left[i]
                tri[1] = i
                tri[2] = right[i]
                triangs.append(tri)
                isPolyPt[i] = False

                # update left and right neighbor lists
                left[right[i]] = left[i]
                right[left[i]] = right[i]

        del right
        del left
        del isPolyPt

        return triangs

    # Original implementation (SKIENA & REVILLA, 2002)
    # creates a ear for the hull
    @staticmethod
    def ear_Q(_i,  _j,  _k,  _p, _isPolyPt):
        t = [None, None, None]  # coordinates for points i,j,k
        t[0] = _p[_i]
        t[1] = _p[_j]
        t[2] = _p[_k]

        # Check for angle ijk (centered in j) greater then 180 degrees
        if Tesselation.cw(t[0], t[1], t[2]):
            return False

        for m in range(0, len(_p)):
            if _isPolyPt[m]:
                if _p[m] != t[0] and _p[m] != t[1] and _p[m] != t[2]:
                    if Tesselation.point_in_triangle(_p[m], t):
                        return False

        return True

    # verifies if the order of the triangle connectivity
    @staticmethod
    def cw(_a,  _b, _c):
        return not (CompGeom.isLeftSide(_a, _b, _c))

    # computes the area of a triangle keeping the sign of orientation
    @staticmethod
    def signed_triangle_area(_a,  _b,  _c):
        return((_a.getX()*_b.getY() - _a.getY()*_b.getX()
                + _a.getY()*_c.getX() - _a.getX()*_c.getY()
                + _b.getX()*_c.getY() - _c.getX()*_b.getY()) / 2.0)

    # verifies if the point _p is inside the triangle _t
    @staticmethod
    def point_in_triangle(_p,  _t):

        for i in range(0, 3):
            if(CompGeom.isRightSide(_t[i], _t[(i+1) % 3], _p)):
                return False
        return True

    @staticmethod
    def tessellate(_pts):
        indices = Tesselation.triangleParing(_pts)
        triangs = []

        for j in range(0, len(indices)):
            triangle = []
            triangle.append(Point(_pts[indices[j][0]].getX(),
                                  _pts[indices[j][0]].getY()))
            triangle.append(Point(_pts[indices[j][1]].getX(),
                                  _pts[indices[j][1]].getY()))
            triangle.append(Point(_pts[indices[j][2]].getX(),
                                  _pts[indices[j][2]].getY()))
            triangs.append(triangle)

        return triangs


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# -------------------------- TOPOLOGICAL ENTITIES ----------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


# Linkedlist superclass declaration
class Linkedlist:
    def __init__(self, prev=None, next=None):

        self.prev = prev
        self.next = next

        if prev is not None:
            self.prev.next = self
        if next is not None:
            self.next.prev = self


# Vertex class declaration
class Vertex(Linkedlist):

    def __init__(self, point=None, he=None):
        Linkedlist.__init__(self)
        self.point = point
        self.he = he
        self.ID = None

    def delete(self):
        if self.next is not None:
            self.next.prev = self.prev
        if self.prev is not None:
            self.prev.next = self.next

    def getType(self):
        return 'VERTEX'

    def incidentFaces(self):
        incFaces = []
        he = self.he
        heBegin = he
        while True:
            face = he.loop.face
            if face not in incFaces:
                incFaces.append(face)

            he = he.mate().next
            if he == heBegin:
                break

        return incFaces

    def incidentEdges(self):
        incEdges = []
        he = self.he
        heBegin = he

        if he.edge is None:
            return incEdges

        while True:
            if he.edge not in incEdges:
                incEdges.append(he.edge)

            he = he.mate().next

            if he == heBegin:
                break

        return incEdges

    def adjacentVertices(self):
        adjVertices = []
        he = self.he
        heBegin = he

        while True:
            he = he.mate()
            if he.vertex != self:
                adjVertices.append(he.vertex)

            he = he.next

            if he == heBegin:
                break

        return adjVertices


class HalfEdge(Linkedlist):

    def __init__(self, vertex=None, loop=None, edge=None, prev=None, next=None):
        Linkedlist.__init__(self, prev, next)
        self.vertex = vertex
        self.edge = edge
        self.loop = loop
        self.ID = None

    # Prepares the half-edge to be deleted
    def delete(self):

        if self.edge is None:
            return
        elif self.next == self:
            self.edge = None
            return self
        else:
            self.edge = None
            self.prev.next = self.next
            self.next.prev = self.prev
            return self.prev

    # Gets the opposite half-edge
    def mate(self):

        if self.edge is None:
            return self.next.prev

        if self == self.edge.he1:
            return self.edge.he2
        else:
            return self.edge.he1

    @staticmethod
    def inBetween(_v1, _v2, _f):

        he = _v1.he
        he_begin = he

        while True:

            if he.mate().vertex == _v2:
                if he.loop.face == _f:
                    return he

            he = he.mate().next

            if he == he_begin:
                return


# Edge class declaration
class Edge(Linkedlist):

    def __init__(self, segment=None, he1=None, he2=None):
        Linkedlist.__init__(self)
        self.he1 = he1
        self.he2 = he2
        self.segment = segment
        self.ID = None

    def delete(self):
        # update linked list
        if self.next is not None:
            self.next.prev = self.prev
        if self.prev is not None:
            self.prev.next = self.next

    def getType(self):
        return 'EDGE'

    def AddHe(self, _v, _where, _sign):

        if _where.edge is None:
            he = _where
        else:
            he = HalfEdge(prev=_where.prev, next=_where)

        he.edge = self
        he.vertex = _v
        he.loop = _where.loop

        if _sign:
            self.he1 = he
        else:
            self.he2 = he

        return he

    def incidentFaces(self):
        incFaces = []
        incFaces.append(self.he1.loop.face)
        if self.he1.loop.face != self.he2.loop.face:
            incFaces.append(self.he2.loop.face)
        return incFaces

    def adjacentEdges(self):
        adjEdges = []
        he1 = self.he1
        he2 = self.he2

        # begin the search at the next he of the first he of the edge
        he = he1.next

        # check if the edge is a closed one
        # if he1 is the half-edge inside the closed edge (do not continue)
        if he != he1:
            while he != he2:
                adjEdges.append(he.edge)
                he = he.mate().next

        # begin the search at the next he of the second he of the edge
        he = he2.next

        # check if the edge is a closed one
        # if he2 is the half-edge inside the closed edge (do not continue)
        if he != he2:
            while he != he1:
                adjEdges.append(he.edge)
                he = he.mate().next

        return adjEdges

    def incidentVertices(self):
        incVertices = []
        incVertices.append(self.he1.vertex)
        incVertices.append(self.he2.vertex)
        return incVertices


# Loop class declaration
class Loop():

    def __init__(self, face=None, he=None, prev=None, next=None):
        self.prev = prev
        self.next = next
        self.face = face
        self.he = he
        self.isClosed = False
        self.ID = None

        if face is not None:
            loopOfFace = self.face.loop
            if loopOfFace is not None:
                self.next = loopOfFace.next
                self.prev = loopOfFace
                loopOfFace.next = self

                if self.next is not None:
                    self.next.prev = self
            else:
                self.face.loop = self

    def delete(self):
        # update linked list
        if self.prev is not None:
            self.prev.next = self.next
        if self.next is not None:
            self.next.prev = self.prev
        if self.face is not None:
            if self == self.face.loop:
                self.face.loop = None


class Face(Linkedlist):

    def __init__(self, shell=None, loop=None, prev=None, next=None, patch=None):
        Linkedlist.__init__(self, prev, next)
        self.shell = shell
        self.loop = loop  # external loop
        self.intLoops = []  # list of internal loops
        self.patch = patch
        self.ID = None

    def delete(self):
        # update linked list
        if self.next is not None:
            self.next.prev = self.prev
        if self.prev is not None:
            self.prev.next = self.next

    def getType(self):
        return 'FACE'

    def adjacentFaces(self):
        adjFaces = []
        loop = self.loop
        if loop.he is not None:
            he = loop.he
            heBegin = he

            while True:
                face = he.mate().loop.face
                if face != self:
                    if face not in adjFaces:
                        adjFaces.append(face)

                he = he.next
                if he == heBegin:
                    break

        return adjFaces

    def incidentEdges(self):
        adjEdges = []
        he = self.loop.he
        heBegin = he

        while True:
            adjEdges.append(he.edge)
            he = he.next

            if he == heBegin:
                break

        return adjEdges

    def incidentVertices(self):
        adjVertexes = []
        he = self.loop.he
        heBegin = he

        while True:
            adjVertexes.append(he.vertex)
            he = he.next

            if he == heBegin:
                break

        return adjVertexes

    def internalFaces(self):
        internalFaces = []

        loop = self.loop.next

        while loop is not None:
            he_begin = loop.he
            he = he_begin

            while True:
                if he.mate().loop != loop:
                    if he.mate().loop.isClosed:
                        internalFaces.append(loop.he.mate().loop.face)
                        break

                he = he.next
                if he == he_begin:
                    break

            loop = loop.next

        return internalFaces

    def updateBoundary(self):
        he_init = self.loop.he
        he = he_init
        bound = []
        orientation = []

        while True:
            if he.edge is not None:
                bound.append(he.edge.segment)
                orientation.append(he == he.edge.he1)

            he = he.next

            if he == he_init:
                break

        self.patch.setBoundary(bound, orientation)

    def updateHoles(self):
        loop = self.loop.next
        self.intLoops.clear()  # clear the list of inner loops
        holes_bound = []
        holes_orientation = []
        intSegments = []
        intSegmentsOrientation = []

        while loop is not None:
            self.intLoops.append(loop)
            he_init = loop.he
            he = he_init
            bound = []
            orientation = []
            isClosed = False

            if he.edge is not None:
                while True:
                    if he.mate().loop.isClosed:
                        isClosed = True

                    bound.append(he.edge.segment)
                    orientation.append(he == he.edge.he1)

                    he = he.next

                    if he == he_init:
                        break

                if len(bound) > 0:
                    if isClosed:
                        holes_bound.append(bound)
                        holes_orientation.append(orientation)
                    else:
                        intSegments.append(bound)
                        intSegmentsOrientation.append(orientation)

            loop = loop.next

        self.patch.setHoles(holes_bound, holes_orientation)
        self.patch.setInternalSegments(intSegments, intSegmentsOrientation)


# Shell class declaration
class Shell:

    def __init__(self, face=None):
        self.face = face
        self.vertices = []
        self.edges = []
        self.faces = []
        self.num_vertices = 0
        self.num_edges = 0
        self.num_faces = -1
        self.num_loops = 0
        self.num_hes = 0

    def insertVertex(self, _vertex):

        if _vertex.ID is None:
            self.num_vertices += 1
            _vertex.ID = self.num_vertices
        elif _vertex.ID > self.num_vertices:
            self.num_vertices = _vertex.ID

        if _vertex.he is not None:
            if _vertex.he.ID is None:
                self.num_hes += 1
                _vertex.he.ID = self.num_hes
            elif _vertex.he.ID > self.num_hes:
                self.num_hes = _vertex.he.ID

            if _vertex.he.loop.ID is None:
                self.num_loops += 1
                _vertex.he.loop.ID = self.num_loops
            elif _vertex.he.loop.ID > self.num_loops:
                self.num_loops = _vertex.he.loop.ID

        if len(self.vertices) > 0:
            _vertex.prev = self.vertices[-1]
            self.vertices[-1].next = _vertex

        self.vertices.append(_vertex)

    def insertEdge(self, _edge):
        if _edge.ID is None:
            self.num_edges += 1
            _edge.ID = self.num_edges
        elif _edge.ID > self.num_edges:
            self.num_edges = _edge.ID

        if _edge.he1 is not None:
            if _edge.he1.ID is None:
                self.num_hes += 1
                _edge.he1.ID = self.num_hes
            elif _edge.he1.ID > self.num_hes:
                self.num_hes = _edge.he1.ID

        if _edge.he2 is not None:
            if _edge.he2.ID is None:
                self.num_hes += 1
                _edge.he2.ID = self.num_hes
            elif _edge.he2.ID > self.num_hes:
                self.num_hes = _edge.he2.ID

        if len(self.edges) > 0:
            _edge.prev = self.edges[-1]
            self.edges[-1].next = _edge

        self.edges.append(_edge)

    def insertFace(self, _face):
        self.faces.append(_face)

        if _face.ID is None:
            self.num_faces += 1
            _face.ID = self.num_faces
        elif _face.ID > self.num_faces:
            self.num_faces = _face.ID

        if _face.loop.ID is None:
            self.num_loops += 1
            _face.loop.ID = self.num_loops
        elif _face.loop.ID > self.num_loops:
            self.num_loops = _face.loop.ID

    def removeVertex(self, _vertex):
        self.vertices.remove(_vertex)

    def removeEdge(self, _edge):
        self.edges.remove(_edge)

    def removeFace(self, _face):
        self.faces.remove(_face)

    def renumberIDS(self):

        self.num_vertices = 0
        self.num_edges = 0
        self.num_faces = -1
        self.num_loops = 0
        self.num_hes = 0

        # renumber vertices IDS
        for vertex in self.vertices:
            self.num_vertices += 1
            vertex.ID = self.num_vertices

        # renumber edges IDS
        for edge in self.edges:
            self.num_edges += 1
            edge.ID = self.num_edges

        # renumber faces,loops and half-edges IDS
        for face in self.faces:
            self.num_faces += 1
            face.ID = self.num_faces
            loop = face.loop

            while loop is not None:
                self.num_loops += 1
                loop.ID = self.num_loops
                he = loop.he
                he_begin = he

                if he is not None:
                    while True:
                        self.num_hes += 1
                        he.ID = self.num_hes

                        he = he.next

                        if he == he_begin:
                            break

                loop = loop.next


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------- EULER OPERATORS ------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


# Note: This code was implemented considering a virtual loop on the infinite face
# (first face to be created) to facilitate the algorithm that migrates the internal
# loops from one face to another. Conceptually this outer loop of the first face does
#  not exist and was just created here to follow the pattern adopted on all faces.
# This virtual loop has no half-edges and therefore no external boundary.


# MakeVertexFaceShell class declaration
class MVFS:
    def __init__(self, point, vertex=None, face=None):

        if point is not None:
            self.shell = Shell()
            self.face = Face(self.shell)
            self.face.patch = Patch()
            self.vertex = Vertex(point)
        else:
            self.vertex = vertex
            self.face = face

    def name(self):
        return 'MVFS'

    def execute(self):
        # create topological entities
        he = HalfEdge(self.vertex)
        loop_out = Loop(self.face)
        new_loop = Loop(self.face, he)

        # set parameters
        self.vertex.he = he
        shell = self.face.shell
        shell.face = self.face
        he.loop = new_loop
        he.prev = he
        he.next = he

    def unexecute(self):
        kvfs = KVFS(self.vertex, self.face)
        kvfs.execute()


# KillVertexFaceShell class declaration
class KVFS:
    def __init__(self, vertex, face):
        self.vertex = vertex
        self.face = face

    def name(self):
        return 'KVFS'

    def execute(self):
        loop_out = self.face.loop
        loop = loop_out.next
        he = self.vertex.he

        self.vertex.he = None
        self.face.loop = None

        self.face.delete()
        self.vertex.delete()
        loop.delete()
        loop_out.delete()
        he.delete()
        del he

    def unexecute(self):
        mvfs = MVFS(None, self.vertex, self.face)
        mvfs.execute()


# MakeEdgeFace class declaration
class MEKR:
    def __init__(self, segment, v_begin, v_end, v_begin_next, v_end_next, face_on, edge=None):

        if segment is not None:
            self.edge = Edge(segment)
        else:
            self.edge = edge

        self.v_begin = v_begin
        self.v_end = v_end
        self.v_begin_next = v_begin_next
        self.v_end_next = v_end_next
        self.face_on = face_on

    def name(self):
        return 'MEKR'

    def execute(self):
        # Loop l2 is always deleted. Therefore, this loop must not be
        # the outter loop of the self.face. This can be avoided by ensuring that
        # he2 is always a half-edge which belongs to an inner loop. This check
        # is currently implemented in previous lines to MEKR usage.
        # [see Class: Hecontroller function: makeEdge for more information.]

        # get the half-edges
        he1 = HalfEdge.inBetween(self.v_begin, self.v_begin_next, self.face_on)
        he2 = HalfEdge.inBetween(self.v_end, self.v_end_next, self.face_on)

        l1 = he1.loop
        l2 = he2.loop

        n_he1 = l2.he

        while True:
            n_he1.loop = l1
            n_he1 = n_he1.next

            if n_he1 == l2.he:
                break

        n_he1 = self.edge.AddHe(he1.vertex, he1, True)
        n_he2 = self.edge.AddHe(he2.vertex, he2, False)
        n_he1.next = he2
        n_he2.next = he1
        he2.prev = n_he1
        he1.prev = n_he2

        l2.delete()

    def unexecute(self):
        kemr = KEMR(self.edge, self.edge.he1.vertex)
        kemr.execute()


# killEdgeFace class declaration
class KEMR:
    def __init__(self, edge, v_out):
        self.edge = edge
        self.v_out = v_out
        self.v_begin = None
        self.v_end = None
        self.v_begin_next = None
        self.v_end_next = None
        self.face_on = None

    def name(self):
        return 'KEMR'

    def execute(self):
        he1 = self.edge.he1
        he2 = self.edge.he2

        if he1.vertex != self.v_out:
            aux = he1
            he1 = he2
            he2 = aux

        # Store the necessary entities for undo
        self.v_begin = he1.vertex
        self.v_end = he2.vertex
        self.face_on = he1.loop.face

        if he2.next == he1:  # Case 1.1 -> he1 vertex is an isolated vertex
            self.v_begin_next = self.v_begin
        elif he2.next.next == he1:  # Case 1.2 -> he1 vertex has a closed segment
            self.v_begin_next = self.v_begin
        else:  # Case 1.3 -> he1 vertex has a branch
            self.v_begin_next = he2.next.mate().vertex

        if he1.next == he2:  # Case 2.1 -> he1 vertex is an isolated vertex
            self.v_end_next = self.v_end
        elif he1.next.next == he2:  # Case 2.2 -> he1 vertex has a closed segment
            self.v_end_next = self.v_end
        else:  # Case 2.3 -> he1 vertex has a branch
            self.v_end_next = he1.next.mate().vertex

        # Now, execute the operation
        ol = he1.loop
        nl = Loop(ol.face)

        he3 = he1.next
        he1.next = he2.next
        he2.next.prev = he1
        he2.next = he3
        he3.prev = he2
        he4 = he2

        while True:
            he4.loop = nl
            he4 = he4.next

            if he4 == he2:
                break

        he1.vertex.he = he1.next
        he2.vertex.he = he2.next

        ol.he = he1.delete()
        nl.he = he2.delete()

        self.edge.he1 = None
        self.edge.he2 = None

        if he1.prev.next != he1:
            del he1
        if he2.prev.next != he2:
            del he2

        self.edge.delete()

    def unexecute(self):
        mekr = MEKR(None, self.v_begin, self.v_end,
                    self.v_begin_next, self.v_end_next, self.face_on, self.edge)
        mekr.execute()


# MakeVertexRing class declaration
class MVR:
    def __init__(self, point, face_on, vertex=None):

        if point is not None:
            self.vertex = Vertex(point)
        else:
            self.vertex = vertex

        self.face_on = face_on

    def name(self):
        return 'MVR'

    def execute(self):
        # create topological entities
        newloop = Loop(self.face_on)
        newhe = HalfEdge(self.vertex, newloop)

        # set parameters
        newloop.he = newhe
        newhe.prev = newhe
        newhe.next = newhe
        self.vertex.he = newhe

    def unexecute(self):
        kvr = KVR(self.vertex, self.face_on)
        kvr.execute()


# KillVertexRing class declaration
class KVR:
    def __init__(self, vertex, face_on):
        self.vertex = vertex
        self.face_on = face_on

    def name(self):
        return 'KVR'

    def execute(self):
        he = self.vertex.he
        loop = he.loop

        self.vertex.he = None

        he.delete()
        loop.delete()

        del he

    def unexecute(self):
        mvr = MVR(None, self.face_on, self.vertex)
        mvr.execute()


# MakeEdgeVertex class declaration
class MEV:
    def __init__(self, point, segment, v_begin, v_next, face_on, vertex=None, edge=None):

        if point is not None:
            self.vertex = Vertex(point)
            self.edge = Edge(segment)
        else:
            self.vertex = vertex
            self.edge = edge

        self.v_begin = v_begin
        self.v_next = v_next
        self.face_on = face_on

    def name(self):
        return 'MEV'

    def execute(self):
        # get half-edges
        he = HalfEdge.inBetween(self.v_begin, self.v_next, self.face_on)

        self.edge.AddHe(he.vertex, he, False)
        self.edge.AddHe(self.vertex, he, True)

        self.vertex.he = he.prev
        he.vertex.he = he

    def unexecute(self):
        kev = KEV(self.edge, self.vertex)
        kev.execute()


# KillEdgeVertex class declaration
class KEV:
    def __init__(self, edge=None, vertex=None):
        self.edge = edge
        self.vertex = vertex
        self.v_begin = None
        self.v_next = None
        self.face_on = None

    def name(self):
        return 'KEV'

    def execute(self):
        he1 = self.edge.he1
        he2 = self.edge.he2

        # switch half-edges such that he1.vertex will be deleted
        if he1.vertex != self.vertex:
            temp = he1
            he1 = he2
            he2 = temp

        # Store the necessary entities for undo
        self.v_begin = he2.vertex

        if he2.next == he1 and he1.next == he2:
            self.v_next = self.v_begin

        elif he2.next != he1 and he1.next == he2:
            self.v_next = he2.next.mate().vertex
        else:
            self.v_next = he1.next.mate().vertex

        self.face_on = he1.loop.face

        # Now, execute the operation
        he = he2.next

        while he != he1:
            he.vertex = he2.vertex
            he = he.mate().next

        he2.vertex.he = he1.next
        he1.loop.he = he1.delete()
        he2.loop.he = he2.delete()

        # cleaning the removed entities
        self.edge.he1 = None
        self.edge.he2 = None
        self.vertex.he = None

        if he1.prev.next != he1:
            del he1

        if he2.prev.next != he2:
            del he2

        self.vertex.delete()
        self.edge.delete()

    def unexecute(self):

        mev = MEV(None, None, self.v_begin, self.v_next,
                  self.face_on, self.vertex, self.edge)
        mev.execute()


# MakeEdgeFace class declaration
class MEF():
    def __init__(self, segment, v_begin, v_end, v_begin_next, v_end_next, face_on, edge=None, face=None):

        if segment is not None:
            self.edge = Edge(segment)
            self.face = Face(face_on.shell)
            self.face.patch = Patch()
        else:
            self.edge = edge
            self.face = face

        self.v_begin = v_begin
        self.v_end = v_end
        self.v_begin_next = v_begin_next
        self.v_end_next = v_end_next
        self.face_on = face_on

    def name(self):
        return 'MEF'

    def execute(self):
        next_face = self.face_on.next
        self.face.prev = self.face_on
        self.face_on.next = self.face
        self.face.next = next_face

        if next_face is not None:
            next_face.prev = self.face

        # get Half-edges
        he1 = HalfEdge.inBetween(
            self.v_begin, self.v_begin_next, self.face_on)
        he2 = HalfEdge.inBetween(
            self.v_end, self.v_end_next, self.face_on)

        newloop = Loop(self.face)
        newloop.isClosed = True

        he = he1
        while he != he2:
            he.loop = newloop
            he = he.next

        nhe1 = self.edge.AddHe(he2.vertex,  he1, False)
        nhe2 = self.edge.AddHe(he1.vertex,  he2, True)

        nhe1.prev.next = nhe2
        nhe2.prev.next = nhe1

        temp = nhe1.prev
        nhe1.prev = nhe2.prev
        nhe2.prev = temp
        newloop.he = nhe1
        nhe1.loop = newloop
        nhe2.loop.he = nhe2

    def unexecute(self):

        kef = KEF(self.edge, self.face)
        kef.execute()


# KillEdgeFace class declaration
class KEF():

    def __init__(self, edge, face):
        self.edge = edge
        self.face = face
        self.v_begin = None
        self.v_end = None
        self.v_begin_next = None
        self.v_end_next = None
        self.face_on = None

    def name(self):
        return 'KEF'

    def execute(self):
        he1 = self.edge.he1
        he2 = self.edge.he2

        # Store the necessary entities for undo
        self.v_begin = he1.vertex
        self.v_end = he2.vertex
        self.v_begin_next = he2.next.next.vertex
        self.v_end_next = he1.next.next.vertex
        self.face_on = he1.loop.face

        # Now, execute the operation
        loop_to_delete = he2.loop
        face_to_delete = self.face
        loop_to_keep = he1.loop

        # set the loops of the half-edges belonging to the loop_to_delete as the loop_to_keep
        he = loop_to_delete.he
        while True:
            he.loop = loop_to_keep
            he = he.next

            if he == loop_to_delete.he:
                break

        # stitch the remaining half-egdes
        he1.prev.next = he2
        he2.prev.next = he1
        he = he2.prev
        he2.prev = he1.prev
        he1.prev = he

        he2.delete()
        he1.delete()

        # set the appropriate half-edges of the vertex that will be kept in the model
        he2.vertex.he = he1.next

        if he2.next != he1:
            he1.vertex.he = he2.next

        # set the appropriate half-edge of the loop_to_keep
        loop_to_keep.he = he1.next

        # reset data of the deleted entities
        self.edge.he1 = None
        self.edge.he2 = None
        self.face.loop = None

        self.edge.delete()
        face_to_delete.delete()
        loop_to_delete.delete()

        if he1.prev.next != he1:
            del he1
        if he2.prev.next != he2:
            del he2

    def unexecute(self):
        mef = MEF(None, self.v_begin, self.v_end, self.v_begin_next,
                  self.v_end_next, self.face_on, self.edge, self.face)
        mef.execute()


# MakeVertexSplitEdge class declaration
class MVSE:

    def __init__(self, point, seg1, seg2, split_edge, vertex=None, edge1=None, edge2=None):

        if point is not None:
            self.vertex = Vertex(point)
            self.edge1 = Edge(seg1)
            self.edge2 = Edge(seg2)
        else:
            self.vertex = vertex
            self.edge1 = edge1
            self.edge2 = edge2

        self.split_edge = split_edge

    def name(self):
        return 'MVSE'

    def execute(self):
        he1 = self.split_edge.he1
        he2 = self.split_edge.he2

        self.edge1.he1 = he1
        self.edge2.he2 = he2
        he1.edge = self.edge1
        he2.edge = self.edge2

        self.edge1.AddHe(self.vertex, he2.next, False)
        self.edge2.AddHe(self.vertex, he1.next, True)
        self.vertex.he = self.edge2.he1

        self.split_edge.he1 = None
        self.split_edge.he2 = None

        self.split_edge.delete()

    def unexecute(self):
        kvje = KVJE(None, self.vertex, self.edge1, self.edge2, self.split_edge)
        kvje.execute()


# KillVertexJoinEdge class declaration
class KVJE:

    def __init__(self, segment, vertex, edge1, edge2, new_edge=None):

        if segment is not None:
            self.new_edge = Edge(segment)
        else:
            self.new_edge = new_edge

        self.vertex = vertex
        self.edge1 = edge1
        self.edge2 = edge2

    def name(self):
        return 'KVJE'

    def execute(self):
        # take the half-edges that are pointing to the vertex
        # that will be deleted
        if self.edge1.he1.vertex == self.vertex:
            he1 = self.edge1.he1
        else:
            he1 = self.edge1.he2

        if self.edge2.he1.vertex == self.vertex:
            he2 = self.edge2.he1
        else:
            he2 = self.edge2.he2

        # adjust the half-edges of the loops because
        # those half-edges will be erased
        he1.loop.he = he1.prev
        he2.loop.he = he2.prev

        old_he1 = he2.prev
        old_he2 = he1.prev

        if he1 == he1.edge.he1:
            self.new_edge.he1 = old_he2
            self.new_edge.he2 = old_he1
        else:
            self.new_edge.he1 = old_he1
            self.new_edge.he2 = old_he2

        he2.delete()
        he1.delete()

        old_he1.edge = self.new_edge
        old_he2.edge = self.new_edge

        self.vertex.he = None
        self.edge1.he1 = None
        self.edge1.he2 = None
        self.edge2.he1 = None
        self.edge2.he2 = None

        self.vertex.delete()
        self.edge1.delete()
        self.edge2.delete()
        del he1
        del he2

    def unexecute(self):
        mvse = MVSE(None, None, None, self.new_edge, self.vertex, self.edge1,
                    self.edge2)
        mvse.execute()


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# -------------------------- AUXILIARY OPERATORS -----------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


# ------------------------------------------------------------

# This code contains a set of operations necessary for
# the proper functioning of the UndoRedo class

# ------------------------------------------------------------

# migrates loops from one face (origin) to a new face (destination)
class MigrateLoops:
    def __init__(self, _origin, _destination, _loops):
        self.origin = _origin
        self.destination = _destination
        self.loops = _loops

    def name(self):
        return 'MIGRATESLOOPS'

    def execute(self):
        for loop_vertex in self.loops:
            loop = self.findLoopOfFace(loop_vertex, self.origin)

            loop.prev.next = loop.next

            if loop.next is not None:
                loop.next.prev = loop.prev

            loop.face = self.destination
            out_loop = self.destination.loop
            loop.prev = out_loop
            loop.next = out_loop.next

            if out_loop.next is not None:
                out_loop.next.prev = loop

            out_loop.next = loop

    def unexecute(self):
        inverse = MigrateLoops(self.destination, self.origin, self.loops)
        inverse.execute()

    def findLoopOfFace(self, _vertex, _face):
        he_begin = _vertex.he
        he = he_begin

        while True:
            if he.loop.face == _face:
                return he.loop

            he = he.mate().next

            if he == he_begin:
                break


# inverts the half-edges of an edge
class Flip:
    def __init__(self, _edge):
        self.edge = _edge

    def name(self):
        return 'FLIP'

    def execute(self):
        temp = self.edge.he1
        self.edge.he1 = self.edge.he2
        self.edge.he2 = temp

    def unexecute(self):
        self.execute()


# This operation turns a face into a hole.
class DelPatch:
    def __init__(self, _patch):
        self.patch = _patch

    def name():
        return 'DEL_PATCH'

    def execute(self):
        self.patch.isDeleted = True
        self.patch.setSelected(False)

    def unexecute(self):
        self.patch.isDeleted = False


# This operation turns a hole into a face.
class CreatePatch:
    def __init__(self, _patch):
        self.patch = _patch

    def name(self):
        return 'CREATE_PATCH'

    def execute(self):
        self.patch.isDeleted = False
        self.patch.setSelected(False)

    def unexecute(self):
        self.patch.isDeleted = True


class InsertShell:
    def __init__(self, _shell, _hemodel):
        self.shell = _shell
        self.hemodel = _hemodel

    def name(self):
        return 'INSERT_SHELL'

    def execute(self):
        self.hemodel.insertShell(self.shell)

    def unexecute(self):
        self.hemodel.removeShell()


class RemoveShell:
    def __init__(self, _shell, _hemodel):
        self.shell = _shell
        self.hemodel = _hemodel

    def name(self):
        return 'REMOVE_SHELL'

    def execute(self):
        self.hemodel.removeShell()

    def unexecute(self):
        self.hemodel.insertShell(self.shell)


class InsertFace:
    def __init__(self, _face, _hemodel):
        self.face = _face
        self.hemodel = _hemodel

    def name(self):
        return 'INSERT_FACE'

    def execute(self):
        self.hemodel.insertFace(self.face)

    def unexecute(self):
        self.hemodel.removeFace(self.face)


class RemoveFace:
    def __init__(self, _face, _hemodel):
        self.face = _face
        self.hemodel = _hemodel

    def name(self):
        return 'REMOVE_FACE'

    def execute(self):
        self.face.patch.setSelected(False)
        self.hemodel.removeFace(self.face)

    def unexecute(self):
        self.hemodel.insertFace(self.face)


class InsertEdge:
    def __init__(self, _edge, _hemodel):
        self.edge = _edge
        self.hemodel = _hemodel

    def name(self):
        return 'INSERT_EDGE'

    def execute(self):
        self.hemodel.insertEdge(self.edge)

    def unexecute(self):
        self.hemodel.removeEdge(self.edge)


class RemoveEdge:
    def __init__(self, _edge, _hemodel):
        self.edge = _edge
        self.hemodel = _hemodel

    def name(self):
        return 'REMOVE_EDGE'

    def execute(self):
        self.edge.segment.setSelected(False)
        self.hemodel.removeEdge(self.edge)

    def unexecute(self):
        self.hemodel.insertEdge(self.edge)


class InsertVertex:
    def __init__(self, _vertex, _hemodel):
        self.vertex = _vertex
        self.hemodel = _hemodel

    def name(self):
        return 'INSERT_VERTEX'

    def execute(self):
        self.hemodel.insertVertex(self.vertex)

    def unexecute(self):
        self.hemodel.removeVertex(self.vertex)


class RemoveVertex:
    def __init__(self, _vertex, _hemodel):
        self.vertex = _vertex
        self.hemodel = _hemodel

    def name(self):
        return 'REMOVE_VERTEX'

    def execute(self):
        self.vertex.point.setSelected(False)
        self.hemodel.removeVertex(self.vertex)

    def unexecute(self):
        self.hemodel.insertVertex(self.vertex)


class SetAttribute:
    def __init__(self, _entity, _attribute):
        self.entity = _entity
        self.attribute = _attribute
        self.oldAttribute = None

        for att in self.entity.attributes:
            if att['type'] == self.attribute['type']:
                self.oldAttribute = att
                break

    def name(self):
        return 'SET_ATTRIBUTE'

    def execute(self):
        if self.oldAttribute is not None:
            self.entity.attributes.remove(self.oldAttribute)
        self.entity.attributes.append(self.attribute)

    def unexecute(self):
        self.entity.attributes.remove(self.attribute)
        if self.oldAttribute is not None:
            self.entity.attributes.append(self.oldAttribute)


class UnSetAttribute:
    def __init__(self, _entity, _attribute):
        self.entity = _entity
        self.attribute = _attribute

    def name(self):
        return 'UNSET_ATTRIBUTE'

    def execute(self):
        self.entity.attributes.remove(self.attribute)

    def unexecute(self):
        self.entity.attributes.append(self.attribute)


class DelAttribute:
    def __init__(self, _attManager, _name, _hemodel):
        self.attManager = _attManager
        self.attribute = self.attManager.getAttributeByName(_name)
        self.hemodel = _hemodel
        self.entities = []

        if self.attribute['applyOnVertex']:
            points = self.hemodel.getPoints()
            for pt in points:
                if self.attribute in pt.attributes:
                    self.entities.append(pt)

        if self.attribute['applyOnEdge']:
            segments = self.hemodel.getSegments()
            for seg in segments:
                if self.attribute in seg.attributes:
                    self.entities.append(seg)

        if self.attribute['applyOnFace']:
            patches = self.hemodel.getPatches()
            for patch in patches:
                if self.attribute in patch.attributes:
                    self.entities.append(patch)

    def name(self):
        return 'DEL_ATTRIBUTE'

    def execute(self):
        self.attManager.removeAttribute(self.attribute)

        for entity in self.entities:
            entity.attributes.remove(self.attribute)

    def unexecute(self):
        attributes = self.attManager.getAttributes()

        for att in attributes:
            if att['name'] == self.attribute['name']:
                self.attribute['name'] = self.attribute['name'] + '_1'

        self.attManager.attributes.append(self.attribute)

        for entity in self.entities:
            setAtt = SetAttribute(entity, self.attribute)
            setAtt.execute()


class SetMesh:
    def __init__(self, _patch, _mesh):
        self.patch = _patch
        self.oldMesh = _patch.mesh
        self.mesh = _mesh

    def name(self):
        return 'SET_MESH'

    def execute(self):
        self.patch.mesh = self.mesh

    def unexecute(self):
        self.patch.mesh = self.oldMesh


class DelMesh:
    def __init__(self, _patch):
        self.patch = _patch
        self.oldMesh = _patch.mesh

    def name(self):
        return 'SET_MESH'

    def execute(self):
        self.patch.mesh = None

    def unexecute(self):
        self.patch.mesh = self.oldMesh


class SetNumberOfSubdivisions:
    def __init__(self, _seg, _attribute):
        self.seg = _seg
        self.numberOfSubdivision = _attribute
        self.oldnumber = _seg.getNumberOfSubdivisions()

    def name(self):
        return 'SET_NUMBER_OF_SUBDIVISIONS'

    def execute(self):
        self.seg.setNumberOfSubdivisions(self.numberOfSubdivision)

    def unexecute(self):
        self.seg.setNumberOfSubdivisions(self.oldnumber)


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ------------------------------ UNDOREDO CLASS ------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


class UndoRedo:

    def __init__(self, limit=-1):
        self.isInserting = False
        self.limit = limit
        self.temp = []
        self.undocommands = []
        self.redocommands = []

    def begin(self):
        if not self.isInserting:
            self.temp = []
            self.isInserting = True

    def end(self):

        if len(self.temp) > 0:

            # insert command
            self.undocommands.insert(0, self.temp)

            # check if reached limit
            if len(self.undocommands) - 1 == self.limit:
                self.undocommands.pop()

            self.clearRedo()

        self.isInserting = False

    def insertOperation(self, _operation):

        if self.isInserting:
            self.temp.insert(0, _operation)

    def lastCommand(self):
        return self.temp

    def lastOperation(self):
        return self.temp[0]

    def hasUndo(self):
        return len(self.undocommands) > 0

    def hasRedo(self):
        return len(self.redocommands) > 0

    def undo(self):
        if not self.isInserting:
            if self.hasUndo():
                self.temp = self.undocommands.pop(0)
                self.redocommands.insert(0, self.temp)

    def redo(self):
        if not self.isInserting:
            if self.hasRedo():
                self.temp = self.redocommands.pop(0)
                self.undocommands.insert(0, self.temp)

    def clear(self):
        self.isInserting = False
        self.clearRedo()
        self.clearUndo()
        self.temp.clear()

    def clearUndo(self):
        self.undocommands.clear()

    def clearRedo(self):
        self.redocommands.clear()


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ------------------------------ HEFILE CLASS --------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


class HeFile():

    @staticmethod
    def saveFile(_shell, _attributes, _filename):

        # get topological entities
        vertices = _shell.vertices
        edges = _shell.edges
        faces = _shell.faces

        # create/ open a file
        split_name = _filename.split('.')
        if split_name[-1] == 'json':
            file = open(f"{_filename}", "w")
        else:
            file = open(f"{_filename}.json", "w")

        # saves the vertices
        vertices_list = []
        for vertex in vertices:

            attributes = vertex.point.attributes
            att_list = []
            for att in attributes:
                att_list.append(att['name'])

            attributes_dict = {
                "att_names": att_list
            }

            if vertex.prev is None:
                prev_ID = None
            else:
                prev_ID = vertex.prev.ID

            if vertex.next is None:
                next_ID = None
            else:
                next_ID = vertex.next.ID

            vertex_dict = {
                'type': 'VERTEX',
                'ID': vertex.ID,
                'prev_ID': prev_ID,
                'next_ID': next_ID,
                'he_ID': vertex.he.ID,
                'point': (vertex.point.getX(), vertex.point.getY()),
                'attributes': attributes_dict
            }
            vertices_list.append(vertex_dict)

        # saves the edges
        edges_list = []
        for edge in edges:

            edge_pts = edge.segment.getPoints()
            pts = []
            for pt in edge_pts:
                pts.append([pt.getX(), pt.getY()])

            attributes = edge.segment.attributes.copy()
            if edge.segment.nsudv is not None:
                attributes.remove(edge.segment.nsudv)
            att_list = []
            for att in attributes:
                att_list.append(att['name'])

            attributes_dict = {
                "nsudv": edge.segment.nsudv,
                "att_names": att_list
            }

            if edge.prev is None:
                prev_ID = None
            else:
                prev_ID = edge.prev.ID

            if edge.next is None:
                next_ID = None
            else:
                next_ID = edge.next.ID

            edge_dict = {
                'type': 'EDGE',
                'ID': edge.ID,
                'prev_ID': prev_ID,
                'next_ID': next_ID,
                'he1_ID': edge.he1.ID,
                'he2_ID': edge.he2.ID,
                'segment_type': f'{edge.segment.getType()}',
                'points': pts,
                'attributes': attributes_dict
            }

            edges_list.append(edge_dict)

        faces_list = []
        for face in faces:

            if face.loop.next is None:
                next_ID = None
            else:
                next_ID = face.loop.next.ID

            # saves the external loop
            he = face.loop.he
            he_begin = he

            if he is None:
                he_list = None
            else:
                he_list = []
                while True:

                    he_dict = {
                        'type': 'HALF-EDGE',
                        'ID': he.ID,
                        'prev_ID': he.prev.ID,
                        'next_ID': he.next.ID,
                        'vertex_ID': he.vertex.ID,
                        'edge_ID': he.edge.ID,
                        'loop_ID': he.loop.ID
                    }

                    he_list.append(he_dict)
                    he = he.next

                    if he == he_begin:
                        break

            loop_dict = {
                'type': 'LOOP',
                'ID': face.loop.ID,
                'prev_ID': None,
                'next_ID': next_ID,
                'face_ID': face.ID,
                'he_loop': he_list,
                'isClosed': face.loop.isClosed

            }

            # saves the internal loops
            intLoops = []
            intLoop = face.loop.next
            while intLoop is not None:

                he = intLoop.he
                he_begin = he

                he_list = []
                while True:

                    if he.edge is None:
                        edge_ID = None
                    else:
                        edge_ID = he.edge.ID

                    he_dict = {
                        'type': 'HALF-EDGE',
                        'ID': he.ID,
                        'prev_ID': he.prev.ID,
                        'next_ID': he.next.ID,
                        'vertex_ID': he.vertex.ID,
                        'edge_ID': edge_ID,
                        'loop_ID': he.loop.ID
                    }

                    he_list.append(he_dict)
                    he = he.next

                    if he == he_begin:
                        break

                if intLoop.next is None:
                    next_ID = None
                else:
                    next_ID = intLoop.next.ID

                intLoop_dict = {
                    'type': 'LOOP',
                    'ID': intLoop.ID,
                    'prev_ID': intLoop.prev.ID,
                    'next_ID': next_ID,
                    'face_ID': face.ID,
                    'he_loop': he_list,
                    'isClosed': intLoop.isClosed
                }

                intLoops.append(intLoop_dict)
                intLoop = intLoop.next

            attributes = face.patch.attributes.copy()
            if face.patch.mesh is not None:
                mesh_dict = face.patch.mesh.mesh_dict
                attributes.remove(mesh_dict)
            else:
                mesh_dict = None

            att_list = []
            for att in attributes:
                att_list.append(att['name'])

            attributes_dict = {
                'isDeleted': face.patch.isDeleted,
                'mesh': mesh_dict,
                "att_names": att_list
            }

            if face.prev is None:
                prev_ID = None
            else:
                prev_ID = face.prev.ID

            if face.next is None:
                next_ID = None
            else:
                next_ID = face.next.ID

            face_dict = {
                'type': 'FACE',
                'ID': face.ID,
                'prev_ID': prev_ID,
                'next_ID': next_ID,
                'loop': loop_dict,
                'intLoops': intLoops,
                'attributes': attributes_dict
            }

            faces_list.append(face_dict)

        shell = {
            'type': 'SHELL',
            'vertices': vertices_list,
            'edges': edges_list,
            'faces': faces_list,
            'attributes_list': _attributes
        }

        json.dump(shell, file, indent=4)
        file.close()

    @ staticmethod
    def loadFile(_file):
        with open(_file, 'r') as file:
            input = json.load(file)

        vertices = input['vertices']
        edges = input['edges']
        faces = input['faces']
        attributes = input['attributes_list']

        # creates the shell
        shell = Shell()

        # creates the edges
        for edge_dict in edges:
            edge = Edge()
            edge.ID = edge_dict['ID']

            # creates a key for the edge
            edge_dict['edge'] = edge

            # set edge segment
            edge_pts = edge_dict['points']
            pts = []
            for pt in edge_pts:
                pts.append(Point(pt[0], pt[1]))

            type = edge_dict['segment_type']

            if type == 'LINE':
                segment = Line(pts[0], pts[1])
            elif type == 'POLYLINE':
                segment = Polyline(pts)

            edge.segment = segment

            # set segment attributes
            att_names = edge_dict['attributes']['att_names']
            for att_name in att_names:
                for attribute in attributes:
                    if att_name == attribute['name']:
                        segment.attributes.append(attribute)

            if edge_dict['attributes']['nsudv'] is not None:
                segment.setNumberOfSubdivisions(
                    edge_dict['attributes']['nsudv'])
                segment.attributes.append(edge_dict['attributes']['nsudv'])

        # creates the vertices
        for vertex_dict in vertices:
            vertex = Vertex()
            vertex.ID = vertex_dict['ID']

            # creates a key for the vertex
            vertex_dict['vertex'] = vertex

            # set the point
            pt = vertex_dict['point']
            vertex.point = Point(pt[0], pt[1])

            # set point attributes
            att_names = vertex_dict['attributes']['att_names']
            for att_name in att_names:
                for attribute in attributes:
                    if att_name == attribute['name']:
                        vertex.point.attributes.append(attribute)

        # creates the faces
        for face_dict in faces:
            face = Face(shell)
            face.patch = Patch()
            face.ID = face_dict['ID']

            # set patch attributes
            att_names = face_dict['attributes']['att_names']
            for att_name in att_names:
                for attribute in attributes:
                    if att_name == attribute['name']:
                        face.patch.attributes.append(attribute)

            # creates a key for the face
            face_dict['face'] = face

            # creates the outer loop
            loop_dict = face_dict['loop']
            loop = Loop(face)
            loop.ID = loop_dict['ID']
            loop.isClosed = loop_dict['isClosed']

            # creates the half-edges
            he_dicts = loop_dict['he_loop']
            if he_dicts is not None:
                for he_dict in he_dicts:
                    he = HalfEdge()
                    he.ID = he_dict['ID']
                    he.loop = loop

                    # creates a key for the he
                    he_dict['he'] = he

                    # set he.vertex and vertex.he
                    for vertex_dict in vertices:
                        if he_dict['vertex_ID'] == vertex_dict['ID']:
                            he.vertex = vertex_dict['vertex']

                            if vertex_dict['he_ID'] == he.ID:
                                he.vertex.he = he

                            break

                    # set he.edge and edge.he(1 or 2)
                    for edge_dict in edges:
                        if he_dict['edge_ID'] == edge_dict['ID']:
                            he.edge = edge_dict['edge']

                            if edge_dict['he1_ID'] == he.ID:
                                he.edge.he1 = he
                                he.edge.segment.setInitPoint(he.vertex.point)
                            else:
                                he.edge.he2 = he
                                he.edge.segment.setEndPoint(he.vertex.point)

                            break

                # set he.prev/next
                he_dicts[0]['he'].prev = he_dicts[-1]['he']
                he_dicts[-1]['he'].next = he_dicts[0]['he']
                for i in range(1, len(he_dicts)):
                    he_dicts[i]['he'].prev = he_dicts[i-1]['he']
                    he_dicts[i-1]['he'].next = he_dicts[i]['he']

                # set loop.he
                loop.he = he_dicts[0]['he']

            # creates internal loops
            intLoops_list = []
            intLoops_dict = face_dict['intLoops']
            for intLoop_dict in intLoops_dict:
                intLoop = Loop()
                intLoop.face = face
                intLoop.ID = intLoop_dict['ID']
                intLoop.isClosed = intLoop_dict['isClosed']
                intLoops_list.append(intLoop)

                # creates the half-edges
                he_dicts = intLoop_dict['he_loop']
                for he_dict in he_dicts:
                    he = HalfEdge()
                    he.ID = he_dict['ID']
                    he.loop = intLoop

                    # creates a key for the he
                    he_dict['he'] = he

                    # set he.vertex and vertex.he
                    for vertex_dict in vertices:
                        if he_dict['vertex_ID'] == vertex_dict['ID']:
                            he.vertex = vertex_dict['vertex']

                            if vertex_dict['he_ID'] == he.ID:
                                he.vertex.he = he

                            break

                    # set he.edge and edge.he(1 or 2)
                    for edge_dict in edges:
                        if he_dict['edge_ID'] == edge_dict['ID']:
                            he.edge = edge_dict['edge']

                            if edge_dict['he1_ID'] == he.ID:
                                he.edge.he1 = he
                                he.edge.segment.setInitPoint(he.vertex.point)
                            else:
                                he.edge.he2 = he
                                he.edge.segment.setEndPoint(he.vertex.point)

                # set he.prev/next
                he_dicts[0]['he'].prev = he_dicts[-1]['he']
                he_dicts[-1]['he'].next = he_dicts[0]['he']
                for i in range(1, len(he_dicts)):
                    he_dicts[i]['he'].prev = he_dicts[i-1]['he']
                    he_dicts[i-1]['he'].next = he_dicts[i]['he']

                # set loop.he
                intLoop.he = he_dicts[0]['he']

            # set loop.prev/next
            if len(intLoops_list) > 0:
                intLoops_list[0].prev = loop
                loop.next = intLoops_list[0]

            for i in range(1, len(intLoops_list)):
                intLoops_list[i].prev = intLoops_list[i-1]
                intLoops_list[i-1].next = intLoops_list[i]

            # set attributes
            attributes_dict = face_dict['attributes']
            face.patch.isDeleted = attributes_dict['isDeleted']

        # set shell face
        shell.face = faces[0]['face']

        # set face prev/next
        for i in range(1, len(faces)):
            faces[i]['face'].prev = faces[i-1]['face']
            faces[i-1]['face'].next = faces[i]['face']

        return vertices, edges, faces, attributes


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ------------------------------ HEMODEL CLASS -------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


class HeModel:

    def __init__(self):
        self.shell = None
        self.infinityFace = None
        self.segments = []
        self.points = []
        self.patches = []
        self.updateSortPatches = False

    def insertShell(self, _shell):
        self.shell = _shell

    def insertVertex(self, _vertex):
        self.shell.insertVertex(_vertex)
        self.points.append(_vertex.point)
        _vertex.point.vertex = _vertex

    def insertEdge(self, _edge):
        self.shell.insertEdge(_edge)
        self.segments.append(_edge.segment)
        _edge.segment.edge = _edge

    def insertFace(self, _face):

        if len(self.shell.faces) == 0:
            self.infinityFace = _face

        self.shell.insertFace(_face)
        _face.patch.face = _face
        self.updateSortPatches = True

    def removeVertex(self, _vertex):
        _vertex.point.vertex = None
        self.shell.removeVertex(_vertex)
        self.points.remove(_vertex.point)

    def removeFace(self, _face):
        if _face == self.infinityFace:
            self.infinityFace = None

        self.shell.removeFace(_face)
        _face.patch.face = None
        self.updateSortPatches = True

    def removeEdge(self, _edge):
        self.shell.removeEdge(_edge)
        self.segments.remove(_edge.segment)
        _edge.segment.edge = None

    def removeShell(self):
        self.shell = None

    def isEmpty(self):
        if self.shell is None:
            return True
        else:
            return False

    def clearAll(self):
        self.shell = None
        self.infinityFace = None
        self.segments = []
        self.points = []
        self.patches = []

    def getPoints(self):
        return self.points

    def getSegments(self):
        return self.segments

    def getPatches(self):

        if self.updateSortPatches:
            self.patches = self.sortPatches()

        return self.patches

    def selectedEdges(self):
        selectedEdges = []

        if self.isEmpty():
            return selectedEdges

        edges = self.shell.edges
        for edge in edges:
            if edge.segment.isSelected():
                selectedEdges.append(edge)

        return selectedEdges

    def selectedVertices(self):

        selectedVertices = []

        if self.isEmpty():
            return selectedVertices

        vertices = self.shell.vertices
        for vertex in vertices:
            if vertex.point.isSelected():
                selectedVertices.append(vertex)

        return selectedVertices

    def selectedFaces(self):

        selectedFaces = []
        if self.isEmpty():
            return selectedFaces

        faces = self.shell.faces
        for face in faces:
            if face.patch.isSelected():
                selectedFaces.append(face)

        return selectedFaces

    def verticesCrossingWindow(self, _xmin, _xmax, _ymin, _ymax):
        vertices = []
        # search the points that are contained in the given rectangle
        vertices_list = self.shell.vertices
        for vertex in vertices_list:
            if _xmin <= vertex.point.getX() and _xmax >= vertex.point.getX():
                if _ymin <= vertex.point.getY() and _ymax >= vertex.point.getY():
                    # then point is in window
                    vertices.append(vertex)

        vertices = list(set(vertices))

        return vertices

    def edgesInWindow(self, _xmin, _xmax, _ymin, _ymax):

        edges_targets = []

        # search the edges that are contained in the given rectangle
        edges_list = self.shell.edges
        for edge in edges_list:
            edge_segment = edge.segment
            edg_xmin, edg_xmax, edg_ymin, edg_ymax = edge_segment.getBoundBox()

            if _xmin <= edg_xmin and _xmax >= edg_xmax:
                if _ymin <= edg_ymin and _ymax >= edg_ymax:
                    # then the edge is in window
                    edges_targets.append(edge)

        return edges_targets

    def edgesCrossingFence(self, _fence):

        edges_targets = []

        xmin, xmax, ymin, ymax = _fence.getBoundBox()

        # get segments crossing fence's bounding box
        edges_list = self.shell.edges
        for edge in edges_list:
            segment = edge.segment
            segment_xmin, segment_xmax, segment_ymin, segment_ymax = segment.getBoundBox()

            if not (xmax < segment_xmin or segment_xmax < xmin or
                    ymax < segment_ymin or segment_ymax < ymin):
                edges_targets.append(edge)

        # Checks if the segment intersects the _fence
        for edge in edges_targets:
            status, pi, param1, param2 = _fence.intersectSegment(edge.segment)

            # If it does not, remove the edge from edgesInFence and go to next edge
            if not status:
                edges_targets.remove(edge)

        return edges_targets

    def edgesCrossingWindow(self, _xmin, _xmax, _ymin, _ymax):
        pts = []

        if _ymin == _ymax or _xmin == _xmax:
            pts.append(Point(_xmin, _ymin))
            pts.append(Point(_xmax, _ymax))
        else:
            # create a retangular fence
            pts.append(Point(_xmin, _ymin))
            pts.append(Point(_xmax, _ymin))
            pts.append(Point(_xmax, _ymax))
            pts.append(Point(_xmin, _ymax))
            pts.append(Point(_xmin, _ymin))

        fence_segment = Polyline(pts)

        edges = self.edgesInWindow(_xmin, _xmax, _ymin, _ymax)

        edges_crossing = self.edgesCrossingFence(fence_segment)
        edges.extend(edges_crossing)

        edges = list(set(edges))  # remove duplicates

        return edges

    def whichFace(self, _pt):
        face = self.infinityFace.next

        while face is not None:
            if face.patch.isPointInside(_pt):
                return face

            face = face.next

        return self.infinityFace

    def sortPatches(self):
        patchesWithoutHoles = []
        facesWithHoles = []

        # initially the faces are organized in two lists of faces with holes
        #  and patches without holes
        faces = self.shell.faces
        for i in range(1, len(faces)):
            if len(faces[i].patch.holes) > 0:
                facesWithHoles.append(faces[i])
            else:
                patchesWithoutHoles.append(faces[i].patch)

        sort_patches = []

        # From this point on, the list of faces with holes is searched looking
        #  for the outermost face. Then the outermost face is added to the new
        #  list of patches with holes
        while len(facesWithHoles) > 0:
            insert = True
            face_target = facesWithHoles[0]
            for j in range(1, len(facesWithHoles)):
                face_point = face_target.loop.he.vertex.point
                poly = facesWithHoles[j].patch.getPoints()

                if CompGeom.isPointInPolygon(poly, face_point):
                    insert = False
                    break

            if insert:
                sort_patches.append(face_target.patch)
                facesWithHoles.pop(0)
            else:
                facesWithHoles.pop(0)
                facesWithHoles.append(face_target)

        sort_patches.extend(patchesWithoutHoles)

        self.updateSortPatches = False

        return sort_patches


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# -------------------------- HECONTROLLER CLASS ------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


class HeController:
    def __init__(self, _hemodel):
        self.undoredo = UndoRedo(10)
        self.attManager = AttribManager()
        self.hemodel = _hemodel
        self.select_segment = True
        self.select_point = True
        self.select_patch = True
        self.file = None
        self.isChanged = False

    def setHeModel(self, _hemodel):
        self.hemodel = _hemodel

    def insertPoint(self, _pt, _tol):
        self.undoredo.begin()

        if type(_pt) == list:
            _pt = Point(_pt[0], _pt[1])

        if self.hemodel.isEmpty():
            shell = self.makeVertexFace(_pt)
            self.hemodel.infinityFace = shell.face
        else:
            self.addPoint(_pt, _tol)

        self.undoredo.end()
        self.update()

    def addPoint(self, _pt, _tol):
        # check whether there is already a point with the same coordinates
        for point in self.hemodel.points:
            tol = Point(_tol, _tol)
            if Point.equal(_pt, point, tol):
                # in this case there is already a vertex with the same coordinates
                return

        # if there isn't one, check whether the point intersects an edge in model
        intersec = False
        edges = self.hemodel.shell.edges
        for edge in edges:
            intersec, param, pi = edge.segment.intersectPoint(_pt, _tol)

            if intersec:
                edge_target = edge
                break

        if intersec:
            # if there is an intersection, then split the edge
            segments = edge_target.segment.splitSegment(param, pi)
            mvse = self.splitEdge(pi, edge_target, segments[0], segments[1])

            # copy edge_split attributes
            mvse.edge1.segment.attributes = edge_target.segment.attributes.copy()
            mvse.edge2.segment.attributes = edge_target.segment.attributes.copy()

        else:
            # if it do not intersect, then find the face where it lies on.
            #  Then add a new vertex to the model
            face_target = self.hemodel.whichFace(_pt)
            self.makeVertexInsideFace(_pt, face_target)

    def insertSegment(self, _segment, _tol):
        self.undoredo.begin()

        if type(_segment) == list:
            pts = []
            while len(_segment) > 0:
                pts.append(Point(_segment.pop(0), _segment.pop(0)))
            _segment = Polyline(pts)

        status, pts, params = _segment.selfIntersect()
        if status:
            # if there are self-intersections, split the segment in segments and
            # then insert each segment at a time
            segment_segments = _segment.split(params, pts)

            for segment in segment_segments:
                if segment is not None:
                    self.addSegment(segment, _tol)
        else:
            self.addSegment(_segment, _tol)

        self.undoredo.end()
        self.update()

    def addSegment(self, _segment, _tol):
        segmentPts = _segment.getPoints()
        init_pt = segmentPts[0]
        end_pt = segmentPts[-1]
        is_closed = (Point.euclidiandistance(init_pt, end_pt) <= _tol)

        if self.hemodel.isEmpty():
            if is_closed:
                # in this case insert the initial point and then the closed segment
                shell = self.makeVertexFace(init_pt)
                self.makeEdge(_segment, init_pt, init_pt)
            else:
                shell = self.makeVertexFace(init_pt)
                self.makeVertexInsideFace(end_pt, shell.face)
                self.makeEdge(_segment, init_pt, end_pt)

        else:
            if is_closed:
                # in this case insert the initial point and then the closed segment
                self.addPoint(init_pt, _tol)

            # intersect incoming edge with existing model
            incoming_edge_split_map, existent_edges_split_map = self.intersectModel(
                _segment, _tol)

            # split the existing edges
            self.splitExistingEdges(existent_edges_split_map)

            # insert incoming segments
            self.insertIncomingSegments(
                _segment, incoming_edge_split_map, _tol)

    def update(self):

        if self.hemodel.isEmpty():
            return

        faces = self.hemodel.shell.faces
        for i in range(1, len(faces)):
            faces[i].updateBoundary()
            faces[i].updateHoles()

        # update internal loops of infinite Face
        loop = faces[0].loop.next
        faces[0].intLoops.clear()
        while loop is not None:
            faces[0].intLoops.append(loop)
            loop = loop.next

        self.isChanged = True

    def makeVertexFace(self, _point):

        # creates, executes and stores the operation
        mvfs = MVFS(_point)
        mvfs.execute()
        self.undoredo.insertOperation(mvfs)

        # insert entities into the model data structure
        insertShell = InsertShell(mvfs.shell, self.hemodel)
        insertShell.execute()
        self.undoredo.insertOperation(insertShell)
        insertFace = InsertFace(mvfs.face, self.hemodel)
        insertFace.execute()
        self.undoredo.insertOperation(insertFace)
        insertVertex = InsertVertex(mvfs.vertex, self.hemodel)
        insertVertex.execute()
        self.undoredo.insertOperation(insertVertex)

        return mvfs

    def makeVertexInsideFace(self, _point, _face):

        # creates, executes and stores the operation
        mvr = MVR(_point, _face)
        mvr.execute()
        self.undoredo.insertOperation(mvr)

        # insert the vertex into the model data structure
        insertVertex = InsertVertex(mvr.vertex, self.hemodel)
        insertVertex.execute()
        self.undoredo.insertOperation(insertVertex)

    def makeEdge(self, _segment, _init_point, _end_point):
        # This function should be used just when the model has already been stitched by the segment.
        # This means that the geometric checks and operations should be done before you call it.
        # If this is the case, then four possibilities may occour:
        # 1: if both end points of the segment belong to vertexes of the model
        # 2: if just the first point of the segment belongs to a vertex of the model
        # 3: if just the last point of the segment belongs to a vertex of the model
        # 4: if the boundary points of the segment do not belong to the model yet

        # check if the points are already present in the model
        initpoint_belongs = False
        endpoint_belongs = False
        init_vertex = _init_point.vertex
        end_vertex = _end_point.vertex

        if init_vertex is not None:
            initpoint_belongs = True

        if end_vertex is not None:
            endpoint_belongs = True

        # update segments points
        _segment.setInitPoint(_init_point)
        _segment.setEndPoint(_end_point)

        if initpoint_belongs and endpoint_belongs:

            begin_tan = _segment.tangent(0.0)
            begin_seg = _segment.curvature(0.0)
            begin_tan = Point.normalize(begin_tan)

            he1 = self.getHalfEdge(init_vertex, begin_tan.getX(), begin_tan.getY(),
                                   -begin_tan.getY(), begin_tan.getX(), begin_seg)

            end_tan = _segment.tangent(1.0)
            end_seg = _segment.curvature(1.0)
            end_tan = Point.normalize(end_tan)

            he2 = self.getHalfEdge(end_vertex, - end_tan.getX(), -end_tan.getY(),
                                   -end_tan.getY(), end_tan.getX(), end_seg)

            if init_vertex.point != end_vertex.point:
                # case 1.1: points are different, then it is an open segment
                # checks if the half-edges have the same loop to decide between MEF and MEKR

                if he1.loop != he2.loop:

                    # case 1.1.1: the half-edges belong to the different loops, then it's a MEKR
                    if he1.loop == he1.loop.face.loop:
                        # if he1 belongs to the outter loop, then no need to inverter
                        mekr = MEKR(_segment, init_vertex, end_vertex, he1.mate(
                        ).vertex, he2.mate().vertex, he1.loop.face)
                        mekr.execute()
                        self.undoredo.insertOperation(mekr)

                    else:
                        # if he2 belongs to the outter loop, then inverter the
                        # half-edges to keep the consistency with the parametric
                        # geometric definition
                        mekr = MEKR(_segment, end_vertex, init_vertex, he2.mate(
                        ).vertex, he1.mate().vertex, he2.loop.face)
                        mekr.execute()
                        self.undoredo.insertOperation(mekr)

                        # inverter the half-edges to keep the consistency with
                        # the parametric geometric definition
                        flip = Flip(mekr.edge)
                        flip.execute()
                        self.undoredo.insertOperation(flip)

                    # insert the entities into the model data structure
                    insertEdge = InsertEdge(mekr.edge, self.hemodel)
                    insertEdge.execute()
                    self.undoredo.insertOperation(insertEdge)

                else:  # case 1.1.2: the half-edges belong to same loops, then it's a MEF

                    existent_loop = he1.loop
                    existent_face = existent_loop.face

                    if self.isSegmentLoopOriented(_segment, he1, he2):
                        mef = MEF(_segment, init_vertex, end_vertex, he1.mate(
                        ).vertex, he2.mate().vertex, existent_face)
                        mef.execute()
                        self.undoredo.insertOperation(mef)
                    else:
                        mef = MEF(_segment, end_vertex, init_vertex, he2.mate(
                        ).vertex, he1.mate().vertex, existent_face)
                        mef.execute()
                        self.undoredo.insertOperation(mef)

                        # inverter the half-edges to keep the consistency with
                        # the parametric geometric definition
                        flip = Flip(mef.edge)
                        flip.execute()
                        self.undoredo.insertOperation(flip)

                    # copy existent_face attributes
                    if existent_face.patch is not None:
                        mef.face.patch.attributes = existent_face.patch.attributes.copy()
                        mef.face.patch.isDeleted = existent_face.patch.isDeleted

                    # insert the entities into the hemodel data structure
                    insertEdge = InsertEdge(mef.edge, self.hemodel)
                    insertEdge.execute()
                    self.undoredo.insertOperation(insertEdge)
                    insertFace = InsertFace(mef.face, self.hemodel)
                    insertFace.execute()
                    self.undoredo.insertOperation(insertFace)

                    mef.face.updateBoundary()

                    inner_loops = self.findInnerLoops(
                        existent_face, mef.face, existent_loop)
                    migrateLoops = MigrateLoops(
                        existent_face, mef.face, inner_loops)
                    migrateLoops.execute()
                    self.undoredo.insertOperation(migrateLoops)

                    mef.face.updateHoles()

            else:
                # case 1.2: points are the same, then it is a closed segment
                split_point = _segment.getPoint(0.5)
                _, param, _ = _segment.intersectPoint(split_point, 0.01)
                seg1, seg2 = _segment.splitSegment(param, split_point)

                if seg1 is None or seg2 is None:
                    print('ERROR: Size of segments are less than geometric tolerance')
                    if len(self.undoredo.temp) > 0:
                        self.undoredo.end()
                        self.undo()
                        self.undoredo.clearRedo()
                    else:
                        self.undoredo.end()
                    raise ValueError

                # --------- Insert first segment ---------------------

                # insert point and segment 1 into the half-edge data structure
                mev = MEV(split_point, seg1, init_vertex, he1.mate().vertex,
                          he1.loop.face)
                mev.execute()
                self.undoredo.insertOperation(mev)

                # # inverter the half-edges to keep the consistency with
                # the parametric geometric definition
                flip = Flip(mev.edge)
                flip.execute()
                self.undoredo.insertOperation(flip)

                insertVertex = InsertVertex(mev.vertex, self.hemodel)
                insertVertex.execute()
                self.undoredo.insertOperation(insertVertex)
                insertEdge1 = InsertEdge(mev.edge, self.hemodel)
                insertEdge1.execute()
                self.undoredo.insertOperation(insertEdge1)

                # --------- Insert second segment ---------------------

                begin_tan = seg2.tangent(0.0)
                begin_seg = seg2.curvature(0.0)
                begin_tan = Point.normalize(begin_tan)

                he1 = self.getHalfEdge(mev.vertex, begin_tan.getX(), begin_tan.getY(),
                                       -begin_tan.getY(), begin_tan.getX(), begin_seg)

                end_tan = seg2.tangent(1.0)
                end_seg = seg2.curvature(1.0)
                end_tan = Point.normalize(end_tan)

                he2 = self.getHalfEdge(end_vertex, - end_tan.getX(), -end_tan.getY(),
                                       -end_tan.getY(), end_tan.getX(), end_seg)

                existent_loop = he1.loop
                existent_face = existent_loop.face

                # check segment orientation
                if self.isSegmentLoopOriented(seg2, he1, he2):
                    mef = MEF(seg2, mev.vertex, end_vertex, he1.mate(
                    ).vertex, he2.mate().vertex, existent_face)
                    mef.execute()
                    self.undoredo.insertOperation(mef)
                else:
                    mef = MEF(seg2, end_vertex, mev.vertex, he2.mate(
                    ).vertex, he1.mate().vertex, existent_face)
                    mef.execute()
                    self.undoredo.insertOperation(mef)

                    # inverter the half-edges to keep the consistency with
                    # the parametric geometric definition
                    flip = Flip(mef.edge)
                    flip.execute()
                    self.undoredo.insertOperation(flip)

                # copy existent_face attributes
                if existent_face.patch is not None:
                    mef.face.patch.attributes = existent_face.patch.attributes.copy()
                    mef.face.patch.isDeleted = existent_face.patch.isDeleted

                # insert the entities into the hemodel data structure
                insertEdge2 = InsertEdge(mef.edge, self.hemodel)
                insertEdge2.execute()
                self.undoredo.insertOperation(insertEdge2)
                insertFace = InsertFace(mef.face, self.hemodel)
                insertFace.execute()
                self.undoredo.insertOperation(insertFace)

                mef.face.updateBoundary()

                inner_loops = self.findInnerLoops(
                    existent_face, mef.face, existent_loop)
                migrateLoops = MigrateLoops(
                    existent_face, mef.face, inner_loops)
                migrateLoops.execute()
                self.undoredo.insertOperation(migrateLoops)

                mef.face.updateHoles()

        elif initpoint_belongs and not endpoint_belongs:
            # case 2: only the initial point of the segment belongs to a vertex of the model

            # get the half-edge of the vertex
            begin_tan = _segment.tangent(0.0)
            begin_seg = _segment.curvature(0.0)
            begin_tan = Point.normalize(begin_tan)
            he = self.getHalfEdge(init_vertex, begin_tan.getX(), begin_tan.getY(),
                                  -begin_tan.getY(), begin_tan.getX(), begin_seg)

            # insert point and incoming segment into the half-edge data structure
            mev = MEV(_end_point, _segment, init_vertex, he.mate().vertex,
                      he.loop.face)
            mev.execute()
            self.undoredo.insertOperation(mev)

            # inverter the half-edges to keep the consistency with
            # the parametric geometric definition
            flip = Flip(mev.edge)
            flip.execute()
            self.undoredo.insertOperation(flip)

            # insert the entities into the model data structure
            insertEdge = InsertEdge(mev.edge, self.hemodel)
            insertEdge.execute()
            self.undoredo.insertOperation(insertEdge)
            insertVertex = InsertVertex(mev.vertex, self.hemodel)
            insertVertex.execute()
            self.undoredo.insertOperation(insertVertex)

        elif not initpoint_belongs and endpoint_belongs:
            # case 3: only the end point of the segment belongs to a vertex of the model

            # get the half-edge of the vertex
            end_tan = _segment.tangent(1.0)
            end_seg = _segment.curvature(1.0)
            end_tan = Point.normalize(end_tan)
            he = self.getHalfEdge(
                end_vertex, - end_tan.getX(), -end_tan.getY(), -end_tan.getY(), end_tan.getX(), end_seg)

            # insert point and incoming segment into the half-edge data structure
            mev = MEV(_init_point, _segment, end_vertex, he.mate().vertex,
                      he.loop.face)
            mev.execute()
            self.undoredo.insertOperation(mev)

            # insert the entities into the model data structure
            insertEdge = InsertEdge(mev.edge, self.hemodel)
            insertEdge.execute()
            self.undoredo.insertOperation(insertEdge)
            insertVertex = InsertVertex(mev.vertex, self.hemodel)
            insertVertex.execute()
            self.undoredo.insertOperation(insertVertex)

        else:
            # case 4: neither of segment's end points belong to the model yet

            # ------------- Insert the init point -------------------

            face_target = self.hemodel.whichFace(_init_point)
            mvr = MVR(_init_point, face_target)
            mvr.execute()
            self.undoredo.insertOperation(mvr)

            # insert the entities into the model data structure
            insertVertex = InsertVertex(mvr.vertex, self.hemodel)
            insertVertex.execute()
            self.undoredo.insertOperation(insertVertex)

            # ----- Insert the point 2 and incoming segment -----

            he = mvr.vertex.he
            mev = MEV(_end_point, _segment, mvr.vertex, he.mate().vertex,
                      he.loop.face)
            mev.execute()
            self.undoredo.insertOperation(mev)

            # inverter the half-edges to keep the consistency with the
            # parametric geometric definition
            flip = Flip(mev.edge)
            flip.execute()
            self.undoredo.insertOperation(flip)

            # insert the entities into the hemodel data structure
            insertEdge = InsertEdge(mev.edge, self.hemodel)
            insertEdge.execute()
            self.undoredo.insertOperation(insertEdge)
            insertVertex = InsertVertex(mev.vertex, self.hemodel)
            insertVertex.execute()
            self.undoredo.insertOperation(insertVertex)

    def delSelectedEntities(self):

        self.undoredo.begin()

        selectedEdges = self.hemodel.selectedEdges()
        selectedVertices = self.hemodel.selectedVertices()

        incidentVertices = []
        for edge in selectedEdges:
            vertices = edge.incidentVertices()
            incidentVertices.extend(vertices)
            self.killEdge(edge)

        incidentVertices = list(set(incidentVertices))  # removes duplicates

        for vertex in incidentVertices:
            if vertex not in selectedVertices:
                self.killVertex(vertex)

        for vertex in selectedVertices:
            edges = vertex.incidentEdges()
            check = False
            if len(edges) == 2:
                check = self.joinEdges(edges[0], edges[1], vertex)

            if not check:
                for edge in edges:
                    vertices = edge.incidentVertices()
                    self.killEdge(edge)

                    for incidentVertex in vertices:
                        if incidentVertex not in selectedVertices:
                            self.killVertex(incidentVertex)

                self.killVertex(vertex)

        selectedFaces = self.hemodel.selectedFaces()
        for face in selectedFaces:
            delPatch = DelPatch(face.patch)
            delPatch.execute()
            self.undoredo.insertOperation(delPatch)

        self.undoredo.end()
        self.update()

    def killVertex(self, _vertex):
        he = _vertex.he

        # # case 1: checks if the vertex which will be deleted belongs
        # to a closed segment (in this case the vertex must not be deleted)
        if he.edge is None:
            # case 1.1 : checks if the vertex which will be deleted is the only one (KVFS)
            vertices = _vertex.he.loop.face.shell.vertices
            if len(vertices) == 1:

                face = _vertex.he.loop.face
                shell = face.shell

                # removes vertex and face from model data structure
                removeFace = RemoveFace(face, self.hemodel)
                removeFace.execute()
                self.undoredo.insertOperation(removeFace)
                removeVertex = RemoveVertex(_vertex, self.hemodel)
                removeVertex.execute()
                self.undoredo.insertOperation(removeVertex)
                removeShell = RemoveShell(shell, self.hemodel)
                removeShell.execute()
                self.undoredo.insertOperation(removeShell)

                # removes shell , face and point from half-edge data structure
                kvfs = KVFS(_vertex, face)
                kvfs.execute()
                self.undoredo.insertOperation(kvfs)

            # case 1.2: the vertex which will be deleted is a floating one (KVR)
            else:

                # removes vertex from model data structure
                removeVertex = RemoveVertex(_vertex, self.hemodel)
                removeVertex.execute()
                self.undoredo.insertOperation(removeVertex)

                kvr = KVR(_vertex, he.loop.face)
                kvr.execute()
                self.undoredo.insertOperation(kvr)

    def killEdge(self, _edge):
        # Case 1: checks if the Edge belongs to a face (its half-edges are
        # in different loops) (KEF)
        # Case 2: checks if both of the Edge's vertexes are incident to more
        #  than one Edge (KEMR)

        he1 = _edge.he1
        he2 = _edge.he2

        if he1.loop != he2.loop:  # Case 1(it's a KEF)

            # find which of its half-edges belongs to an outter loop
            if he1.loop == he1.loop.face.loop:
                face_to_delete = he1.loop.face
                face_to_keep = he2.loop.face
            else:
                face_to_delete = he2.loop.face
                face_to_keep = he1.loop.face

            # store inner loops
            loop = face_to_delete.loop.next
            inner_loops = []

            while loop is not None:
                inner_loops.append(loop.he.vertex)
                loop = loop.next

            # migrate loops to face_to_keep
            migrateLoops = MigrateLoops(
                face_to_delete, face_to_keep, inner_loops)
            migrateLoops.execute()
            self.undoredo.insertOperation(migrateLoops)

            # checks if it is necessary to invert the half-edges
            if he1.loop.face == face_to_delete:
                # inverter the half-edges to keep the consistency with
                # the parametric geometric definition
                flip = Flip(_edge)
                flip.execute()
                self.undoredo.insertOperation(flip)

            # removes face_to_delete from model data structure
            removeFace = RemoveFace(face_to_delete, self.hemodel)
            removeFace.execute()
            self.undoredo.insertOperation(removeFace)
            removeEdge = RemoveEdge(_edge, self.hemodel)
            removeEdge.execute()
            self.undoredo.insertOperation(removeEdge)

            kef = KEF(_edge, face_to_delete)
            kef.execute()
            self.undoredo.insertOperation(kef)

        else:
            # Test whether the edge belongs to the outter loop of its face
            vertex_out = he1.vertex
            if he1.loop == he1.loop.face.loop:
                if self.isLoopCCW(he1.next, he2):
                    vertex_out = he2.vertex

                    # inverter the half-edges to keep the consistency with
                    # the parametric geometric definition
                    flip = Flip(_edge)
                    flip.execute()
                    self.undoredo.insertOperation(flip)

            # removes edge from model data structure
            removeEdge = RemoveEdge(_edge, self.hemodel)
            removeEdge.execute()
            self.undoredo.insertOperation(removeEdge)

            # removes edge from half-edge data structure
            kemr = KEMR(_edge, vertex_out)
            kemr.execute()
            self.undoredo.insertOperation(kemr)

    def getHalfEdge(self, vertex, _tanx, _tany, _normx, _normy, _curvature):

        # get the incident edges of the vertex
        edges = vertex.incidentEdges()

        # case the vertex contains only one edge then returns its half-edge,
        # otherwise returns the half-edge that is most right of the "new edge"
        if len(edges) < 2:
            return vertex.he

        # computes the angle with the horizontal for the "new edge"
        angle_min = 2*CompGeom.PI
        curv_vec_norm_min = 0
        curv_vec_norm_i = 0
        curv_vec_norm_min_first = True
        angleRef = math.atan2(_tany, _tanx)

        if angleRef < 0:
            angleRef += 2*CompGeom.PI

        # find vector normal to given tangent
        ref_norm = Point.normalize(Point(-_tany, _tany))
        curv_vec_ref = Point(_normx*_curvature, _normy*_curvature)
        dotprod_ref = Point.dotprod(curv_vec_ref, ref_norm)

        # loops over the vertex edges to identify the desired half-edge
        he_i = vertex.he

        while True:
            # computes the angle with the horizontal for the "current edge"
            # get the correct tangent

            if he_i == he_i.edge.he1:
                tan = Point.normalize(he_i.edge.segment.tangent(0.0))
                segment_curvature = he_i.edge.segment.curvature(0.0)
                curv_vec_i = Point(-tan.getY() * segment_curvature,
                                   tan.getX() * segment_curvature)
                angle_i = math.atan2(tan.getY(), tan.getX())
            else:
                tan = Point.normalize(he_i.edge.segment.tangent(1.0))
                segment_curvature = he_i.edge.segment.curvature(1.0)
                curv_vec_i = Point(-tan.getY() * segment_curvature,
                                   tan.getX() * segment_curvature)
                angle_i = math.atan2(-tan.getY(), -tan.getX())

            if angle_i < 0:
                angle_i += 2 * CompGeom.PI

            # obtains only positive values from reference edge in ccw
            angle_i = angleRef - angle_i

            if angle_i < 0:
                angle_i = angle_i + 2.0 * CompGeom.PI

            # check if model segment is above incoming
            if angle_i == 0.0 and Point.dotprod(curv_vec_i, ref_norm) > dotprod_ref:
                angle_i = 2.0 * CompGeom.PI

            if angle_i < angle_min:
                angle_min = angle_i
                he_min = he_i
            elif angle_i == angle_min:  # tie break using curvature
                curv_vec_norm_i = Point.dotprod(curv_vec_i, curv_vec_i)

                if curv_vec_norm_min_first:
                    curv_vec_norm_min_first = False
                    curv_vec_norm_min = curv_vec_norm_i
                elif curv_vec_norm_i < curv_vec_norm_min:
                    curv_vec_norm_min = curv_vec_norm_i
                    he_min = he_i

            he_i = he_i.mate().next

            if he_i == vertex.he:
                break

        return he_min

    def intersectModel(self, _segment, _tol):
        incoming_edge_split_map = []
        existent_edges_split_map = []

        # gets the incoming segment bounding box
        xmin, xmax, ymin, ymax = _segment.getBoundBox()

        # -------------------------VERTEX INTERSECTION-------------------------
        # OBS: only floating vertices
        verticesInBound = self.hemodel.verticesCrossingWindow(
            xmin, xmax, ymin, ymax)
        for vertex in verticesInBound:
            if vertex.he.edge is None:
                status, param, pi = _segment.intersectPoint(
                    vertex.point, _tol)
                if status:
                    incoming_edge_split_map.append([param, vertex.point])

        # -------------------------EDGE INTERSECTION---------------------------
        edgesInBound = self.hemodel.edgesCrossingWindow(
            xmin, xmax, ymin, ymax)
        for edge in edgesInBound:
            existent_edge_split_map = []
            segment = edge.segment
            status, pts, existent_params, incoming_params = segment.intersectSegment(
                _segment)

            if status:
                for i in range(0, len(pts)):
                    if abs(existent_params[i]) <= CompGeom.ABSTOL:
                        point = edge.he1.vertex.point
                    elif abs(existent_params[i]-1.0) <= CompGeom.ABSTOL:
                        point = edge.he2.vertex.point
                    else:
                        point = pts[i]
                        # insert at existent params map
                        existent_edge_split_map.append(
                            [existent_params[i], point])

                    # insert in incoming params map
                    incoming_edge_split_map.append(
                        [incoming_params[i], point])

                if len(existent_edge_split_map) > 0:

                    # removes duplicate elements
                    uniqueList = []
                    for item in existent_edge_split_map:
                        insert = True
                        for unique_item in uniqueList:
                            if abs(item[0]-unique_item[0]) <= _tol:
                                tol = Point(_tol, _tol)
                                if Point.equal(item[1], unique_item[1], tol):
                                    insert = False
                                    break

                        if insert:
                            uniqueList.append(item)

                    existent_edge_split_map = uniqueList
                    existent_edge_split_map.sort()

                    existent_edges_split_map.append(
                        [edge, existent_edge_split_map])

        # removes duplicate elements
        uniqueList = []
        for item in incoming_edge_split_map:
            if item not in uniqueList:
                uniqueList.append(item)

        incoming_edge_split_map = uniqueList
        incoming_edge_split_map.sort()

        # try to insert init and end points
        segment_pts = _segment.getPoints()
        if len(incoming_edge_split_map) == 0:
            incoming_edge_split_map.append([0.0, segment_pts[0]])
            incoming_edge_split_map.append([1.0, segment_pts[-1]])
        else:
            if incoming_edge_split_map[0][0] != 0.0:
                incoming_edge_split_map.insert(0, [0.0, segment_pts[0]])
            if incoming_edge_split_map[-1][0] != 1.0:
                incoming_edge_split_map.append([1.0, segment_pts[-1]])

        return incoming_edge_split_map, existent_edges_split_map

    def splitExistingEdges(self, _edges_split_map):

        # split each intersected existent segment and insert its segments
        for edge_split_map in _edges_split_map:
            # geometrically split segments
            split_params = []
            split_pts = []
            existent_edge = edge_split_map[0]
            for split_nodes in edge_split_map[1]:
                split_params.append(split_nodes[0])
                split_pts.append(split_nodes[1])

            segments = existent_edge.segment.split(split_params, split_pts)

            # for each split point, split the segment and insert seg1
            # seg2 will have the correct topology of the splitted existent edge
            # and the geometric information of the splitted existent edge
            # each subsequent call will insert a segment (seg1) which will both
            # have geometric and topological information
            # this loop go as far as there is more than 2 segments remaining

            initial_segment = existent_edge.segment.clone()
            while len(segments) > 2:

                # split the existent segment
                segment1, segment2 = initial_segment.splitSegment(
                    split_params[0], split_pts[0])

                # split the edge
                mvse = self.splitEdge(split_pts[0],
                                      existent_edge, segments[0], segment2)

                # copy edge_split attributes
                mvse.edge1.segment.attributes = existent_edge.segment.attributes.copy()
                mvse.edge2.segment.attributes = existent_edge.segment.attributes.copy()

                # update the next segment to be splitted
                existent_edge = mvse.edge2

                segments.pop(0)
                split_params.pop(0)
                split_pts.pop(0)

            # at this point there are only two segments to be inserted
            # then insert them both at the same time
            mvse = self.splitEdge(split_pts[0], existent_edge,
                                  segments[0], segments[1])

            # copy edge_split attributes
            mvse.edge1.segment.attributes = existent_edge.segment.attributes.copy()
            mvse.edge2.segment.attributes = existent_edge.segment.attributes.copy()

    def splitEdge(self, _pt, _split_edge, _seg1, _seg2):

        if _seg1 is None or _seg2 is None:
            if len(self.undoredo.temp) > 0:
                self.undoredo.end()
                self.undo()
                self.undoredo.clearRedo()
            else:
                self.undoredo.end()
            print('Unable to insert the segment')
            raise ValueError

        # insert and remove entities in model data structure
        removeEdge = RemoveEdge(_split_edge, self.hemodel)
        removeEdge.execute()
        self.undoredo.insertOperation(removeEdge)

        mvse = MVSE(_pt, _seg1, _seg2, _split_edge)
        mvse.execute()
        self.undoredo.insertOperation(mvse)

        insertVertex = InsertVertex(mvse.vertex, self.hemodel)
        insertVertex.execute()
        self.undoredo.insertOperation(insertVertex)
        insertEdge = InsertEdge(mvse.edge1, self.hemodel)
        insertEdge.execute()
        self.undoredo.insertOperation(insertEdge)
        insertEdge = InsertEdge(mvse.edge2, self.hemodel)
        insertEdge.execute()
        self.undoredo.insertOperation(insertEdge)

        return mvse

    def joinEdges(self, _edge1, _edge2, _vertex):

        loop1 = _edge1.he1.loop
        loop2 = _edge1.he2.loop

        if self.checkClosedSegment(loop1) or self.checkClosedSegment(loop2):
            return False

        seg1_pts = _edge1.segment.getPoints().copy()
        seg2_pts = _edge2.segment.getPoints().copy()
        joinned_pts = []

        if seg1_pts[0] == _vertex.point:
            init_pt1 = True
        else:
            init_pt1 = False

        if seg2_pts[0] == _vertex.point:
            init_pt2 = True
        else:
            init_pt2 = False

        if init_pt1 and init_pt2:
            seg1_pts.reverse()
            seg1_pts.pop()
            joinned_pts.extend(seg1_pts)
            joinned_pts.extend(seg2_pts)
        elif not init_pt1 and not init_pt2:
            seg1_pts.pop()
            joinned_pts.extend(seg1_pts)
            seg2_pts.reverse()
            joinned_pts.extend(seg2_pts)
        elif init_pt1 and not init_pt2:
            joinned_pts.extend(seg2_pts)
            joinned_pts.pop()
            joinned_pts.extend(seg1_pts)
        elif not init_pt1 and init_pt2:
            joinned_pts.extend(seg1_pts)
            joinned_pts.pop()
            joinned_pts.extend(seg2_pts)

        segment = Polyline(joinned_pts)

        # removes entities
        removeVertex = RemoveVertex(_vertex, self.hemodel)
        removeVertex.execute()
        self.undoredo.insertOperation(removeVertex)

        removeEdge = RemoveEdge(_edge1, self.hemodel)
        removeEdge.execute()
        self.undoredo.insertOperation(removeEdge)

        removeEdge = RemoveEdge(_edge2, self.hemodel)
        removeEdge.execute()
        self.undoredo.insertOperation(removeEdge)

        # execute Euler operation
        if _edge1.he1.vertex == _vertex:
            if _edge2.he1.vertex == _vertex:
                flip = Flip(_edge1)
                flip.execute()
                self.undoredo.insertOperation(flip)

                kvje = KVJE(segment, _vertex, _edge1, _edge2)
                kvje.execute()
                self.undoredo.insertOperation(kvje)
            else:
                kvje = KVJE(segment, _vertex, _edge2, _edge1)
                kvje.execute()
                self.undoredo.insertOperation(kvje)
        else:
            if not _edge2.he1.vertex == _vertex:
                flip = Flip(_edge2)
                flip.execute()
                self.undoredo.insertOperation(flip)

            kvje = KVJE(segment, _vertex, _edge1, _edge2)
            kvje.execute()
            self.undoredo.insertOperation(kvje)

        # copy attributes
        if len(kvje.edge1.segment.attributes) > 0:
            kvje.new_edge.segment.attributes = kvje.edge1.segment.attributes.copy()
        else:
            kvje.new_edge.segment.attributes = kvje.edge2.segment.attributes.copy()

        # insert joinned edge in model data structure
        insertEdge = InsertEdge(kvje.new_edge, self.hemodel)
        insertEdge.execute()
        self.undoredo.insertOperation(insertEdge)

        return True

    def checkClosedSegment(self, _loop):
        he_begin = _loop.he
        he = he_begin

        vertices = []
        while True:
            if he.vertex not in vertices:
                vertices.append(he.vertex)
            he = he.next

            if he == he_begin:
                break

        if len(vertices) == 2:
            return True
        else:
            return False

    def insertIncomingSegments(self, _segment, _incoming_segment_split_map, _tol):
        # get the splitted segments
        split_params = []
        split_pts = []
        points = []
        tol = Point(_tol, _tol)

        for split_nodes in _incoming_segment_split_map:
            split_params.append(split_nodes[0])
            split_pts.append(split_nodes[1])
            points.append(split_nodes[1])

        split_params.pop(0)
        split_params.pop()
        split_pts.pop(0)
        split_pts.pop()

        segments = _segment.split(split_params, split_pts)

        # insert segments in model
        init_point = points.pop(0)

        for seg in segments:

            if seg is None:
                if len(self.undoredo.temp) > 0:
                    self.undoredo.end()
                    self.undo()
                    self.undoredo.clearRedo()
                else:
                    self.undoredo.end()
                print('it was not possible to insert the segment')
                raise ValueError

            # get end vertex and increment
            end_point = points.pop(0)

            # The list of vertices of the hemodel is checked, verifying if the
            # init_point and end_point are already exists in the model
            init_vertex = None
            end_vertex = None
            vertices = self.hemodel.shell.vertices
            for vertex in vertices:
                if Point.equal(vertex.point, init_point, tol):
                    init_vertex = vertex
                    init_point = init_vertex.point

                if Point.equal(vertex.point, end_point, tol):
                    end_vertex = vertex
                    end_point = end_vertex.point

            make_segment = True
            if seg.length(0, 1) <= _tol:
                make_segment = False
            # check if the segment to be inserted already exists in the model
            elif init_vertex is not None and end_vertex is not None:
                if init_vertex.he is not None and end_vertex.he is not None:
                    edgesBetween = self.edgesBetween(init_vertex, end_vertex)
                    for edge in edgesBetween:
                        if seg.isEqual(edge.segment, _tol):
                            make_segment = False
                            break

            # insert the segment
            if make_segment:
                self.makeEdge(seg, init_point, end_point)

            # change the initial point
            init_point = end_point

    def isSegmentLoopOriented(self, _segment, _he1, _he2):
        he = _he1
        area = 0.0
        cont = 0
        while he != _he2:
            if he == he.edge.he1:
                area += he.edge.segment.boundIntegral()
            else:
                area -= he.edge.segment.boundIntegral()

            he = he.next

        area -= _segment.boundIntegral()

        return area >= 0

    def isLoopCCW(self, _he1, _he2):
        area = 0.0
        he = _he1

        while he != _he2:
            if he == he.edge.he1:
                area += he.edge.segment.boundIntegral()
            else:
                area -= he.edge.segment.boundIntegral()

            he = he.next

        return area > CompGeom.ABSTOL

    def findInnerLoops(self, _existent_face, _new_face, _existent_loop):
        loop = _existent_face.loop.next
        inner_loops = []

        while loop is not None:
            if loop != _existent_loop:
                if _new_face.patch.isPointInside(loop.he.vertex.point):
                    inner_loops.append(loop.he.vertex)

            loop = loop.next

        return inner_loops

    def edgesBetween(self, _v1, _v2):
        segments_between = []
        he = _v1.he
        he_begin = he

        # check for floating vertex
        if he.edge is None:
            return segments_between

        while True:
            he = he.mate()
            if he.vertex == _v2:
                segments_between.append(he.edge)

            he = he.next

            if he == he_begin:
                break

        return segments_between

    def createPatch(self):
        self.undoredo.begin()
        selectedFaces = self.hemodel.selectedFaces()

        for face in selectedFaces:
            if face.patch.isDeleted:
                createPatch = CreatePatch(face.patch)
                createPatch.execute()
                self.undoredo.insertOperation(createPatch)

        self.isChanged = True
        self.undoredo.end()

    def undo(self):
        # check whether has redo
        if not self.undoredo.hasUndo():
            return

        # update undo stack
        self.undoredo.undo()

        # undo last command
        lastCommand = self.undoredo.lastCommand()
        for comand in lastCommand:
            comand.unexecute()

        self.update()

    def redo(self):
        # check whether has redo
        if not self.undoredo.hasRedo():
            return

        # update redo stack
        self.undoredo.redo()

        lastCommand = self.undoredo.lastCommand()
        for i in range(len(lastCommand)-1, -1, -1):
            lastCommand[i].execute()

        self.update()

    def selectPick(self, _x,  _y,  _tol,  _shiftkey):

        if self.hemodel.isEmpty():
            return

        # select point
        ispointSelected = False
        id_target = -1
        dmin = _tol
        points = self.hemodel.getPoints()
        if self.select_point:
            for i in range(0, len(points)):
                dist = Point.euclidiandistance(Point(_x, _y), points[i])
                if dist < dmin:
                    dmin = dist
                    id_target = i

            # Revert selection of picked point
            if id_target > -1:
                ispointSelected = True
                if points[id_target].isSelected():
                    points[id_target].setSelected(False)
                else:
                    points[id_target].setSelected(True)

        if not _shiftkey:
            # If shift key is not pressed, unselect all points except
            # the picked one (if there was one selected)
            for i in range(0, len(points)):
                if i != id_target:
                    points[i].setSelected(False)

        # select segment
        issegmentselected = False
        id_target = -1
        dmin = _tol
        segments = self.hemodel.getSegments()
        if self.select_segment and not ispointSelected:
            for i in range(0, len(segments)):
                # Compute distance between given point and segment and
                # update minimum distance
                xC, yC, d = segments[i].closestPoint(_x, _y)
                if d < dmin:
                    dmin = d
                    id_target = i

            # Revert selection of picked segment
            if id_target > -1:
                issegmentselected = True
                if segments[id_target].isSelected():
                    segments[id_target].setSelected(False)
                else:
                    segments[id_target].setSelected(True)

        if not _shiftkey:
            # If shift key is not pressed, unselect all segments except
            # the picked one (if there was one selected)
            for i in range(0, len(segments)):
                if i != id_target:
                    segments[i].setSelected(False)

        patches = self.hemodel.getPatches()
        if self.select_patch and not ispointSelected and not issegmentselected:
            # Check whether point is inside a patch
            p = Point(_x, _y)

            for i in range(0, len(patches)):
                if patches[i].isPointInside(p):
                    if patches[i].isSelected():
                        patches[i].setSelected(False)
                    else:
                        patches[i].setSelected(True)
                else:
                    if not _shiftkey:
                        patches[i].setSelected(False)
        elif not _shiftkey:
            for i in range(0, len(patches)):
                patches[i].setSelected(False)

    def selectFence(self, _xmin, _xmax, _ymin, _ymax, _shiftkey):

        if self.hemodel.isEmpty():
            return

        segments = self.hemodel.getSegments()
        if self.select_segment:
            # select segments
            for i in range(0, len(segments)):
                xmin_c, xmax_c, ymin_c, ymax_c = segments[i].getBoundBox(
                )
                if ((xmin_c < _xmin) or (xmax_c > _xmax) or
                        (ymin_c < _ymin) or (ymax_c > _ymax)):
                    inFence = False
                else:
                    inFence = True

                if inFence:
                    # Select segment inside fence
                    segments[i].setSelected(True)
                else:
                    if not _shiftkey:
                        segments[i].setSelected(False)
        elif not _shiftkey:
            for i in range(0, len(segments)):
                segments[i].setSelected(False)

        points = self.hemodel.getPoints()
        if self.select_point:
            # select points
            for i in range(0, len(points)):
                x = points[i].getX()
                y = points[i].getY()

                if ((x < _xmin) or (x > _xmax) or
                        (y < _ymin) or (y > _ymax)):
                    inFence = False
                else:
                    inFence = True

                if inFence:
                    # Select segment inside fence
                    points[i].setSelected(True)
                else:
                    if not _shiftkey:
                        points[i].setSelected(False)
        elif not _shiftkey:
            for i in range(0, len(points)):
                points[i].setSelected(False)

        patches = self.hemodel.getPatches()
        if self.select_patch:
            # select patches
            for i in range(0, len(patches)):
                xmin_r, xmax_r, ymin_r, ymax_r = patches[i].getBoundBox(
                )
                if((xmin_r < _xmin) or (xmax_r > _xmax) or
                        (ymin_r < _ymin) or (ymax_r > _ymax)):
                    inFence = False
                else:
                    inFence = True

                if inFence:
                    # Select patch inside fence
                    patches[i].setSelected(True)
                else:
                    # If shift key is not pressed, unselect patch outside fence
                    if not _shiftkey:
                        patches[i].setSelected(False)
        elif not _shiftkey:
            for i in range(0, len(patches)):
                patches[i].setSelected(False)

    def unSelectAll(self):
        points = self.hemodel.getPoints()
        segments = self.hemodel.getSegments()
        patches = self.hemodel.getPatches()

        for point in points:
            point.setSelected(False)

        for segment in segments:
            segment.setSelected(False)

        for patch in patches:
            patch.setSelected(False)

    def changePointSelect(self, _select):
        self.select_point = _select

    def changeSegmentSelect(self, _select):
        self.select_segment = _select

    def changePatchSelect(self, _select):
        self.select_patch = _select

    def saveFile(self, _filename):

        self.file = _filename

        # renumber IDS
        self.hemodel.shell.renumberIDS()

        shell = self.hemodel.shell
        HeFile.saveFile(shell, self.attManager.getAttributes(), _filename)
        self.isChanged = False

    def openFile(self, _filename):

        self.file = _filename
        vertices, edges, faces, attributes = HeFile.loadFile(_filename)

        self.undoredo.clear()
        self.hemodel.clearAll()

        shell = faces[0]['face'].shell
        self.hemodel.insertShell(shell)
        self.attManager.attributes = attributes

        for vertex_dict in vertices:
            self.hemodel.insertVertex(vertex_dict['vertex'])

        for edge_dict in edges:
            self.hemodel.insertEdge(edge_dict['edge'])

        for face_dict in faces:
            self.hemodel.insertFace(face_dict['face'])

        self.update()
        self.isChanged = False

    def setAttribute(self, _name):

        attribute = self.attManager.getAttributeByName(_name)
        self.undoredo.begin()

        if attribute['applyOnVertex']:
            points = self.hemodel.getPoints()

            for pt in points:
                if pt.isSelected():
                    pt.setSelected(False)
                    setAtt = SetAttribute(pt, attribute)
                    setAtt.execute()
                    self.undoredo.insertOperation(setAtt)

        if attribute['applyOnEdge']:
            segments = self.hemodel.getSegments()

            for seg in segments:
                if seg.isSelected():
                    seg.setSelected(False)
                    setAtt = SetAttribute(seg, attribute)
                    setAtt.execute()
                    self.undoredo.insertOperation(setAtt)

                    # change the support conditions of the segment points
                    if attribute['type'] == 'Support Conditions':
                        seg_vertices = seg.edge.incidentVertices()
                        if attribute not in seg_vertices[0].point.attributes:
                            setAtt = SetAttribute(
                                seg_vertices[0].point, attribute)
                            setAtt.execute()
                            self.undoredo.insertOperation(setAtt)

                        if attribute not in seg_vertices[1].point.attributes:
                            setAtt = SetAttribute(
                                seg_vertices[-1].point, attribute)
                            setAtt.execute()
                            self.undoredo.insertOperation(setAtt)

        if attribute['applyOnFace']:
            patches = self.hemodel.getPatches()

            for patch in patches:
                if patch.isSelected() and not patch.isDeleted:
                    patch.setSelected(False)
                    setAtt = SetAttribute(patch, attribute)
                    setAtt.execute()
                    self.undoredo.insertOperation(setAtt)

        self.undoredo.end()
        self.isChanged = True

    def unSetAttribute(self, _name):

        attribute = self.attManager.getAttributeByName(_name)
        self.undoredo.begin()

        if attribute['applyOnVertex']:
            points = self.hemodel.getPoints()

            for pt in points:
                if pt.isSelected():
                    pt.setSelected(False)
                    if attribute in pt.attributes:
                        unsetAtt = UnSetAttribute(pt, attribute)
                        unsetAtt.execute()
                        self.undoredo.insertOperation(unsetAtt)

        if attribute['applyOnEdge']:
            segments = self.hemodel.getSegments()

            for seg in segments:
                if seg.isSelected():
                    seg.setSelected(False)
                    if attribute in seg.attributes:
                        unsetAtt = UnSetAttribute(seg, attribute)
                        unsetAtt.execute()
                        self.undoredo.insertOperation(unsetAtt)

                        # update mesh
                        if attribute['type'] == 'Number of Subdivisions':
                            face1 = seg.edge.he1.loop.face
                            face2 = seg.edge.he2.loop.face

                            if face1.patch.mesh is not None:
                                self.delMesh(face1)

                            if face2.patch.mesh is not None:
                                self.delMesh(face2)

                        # change the support conditions of the segment points
                        elif attribute['type'] == 'Support Conditions':
                            seg_vertices = seg.edge.incidentVertices()
                            if attribute in seg_vertices[0].point.attributes:
                                unsetAtt = UnSetAttribute(
                                    seg_vertices[0].point, attribute)
                                unsetAtt.execute()
                                self.undoredo.insertOperation(unsetAtt)

                            if attribute in seg_vertices[1].point.attributes:
                                unsetAtt = UnSetAttribute(
                                    seg_vertices[-1].point, attribute)
                                unsetAtt.execute()
                                self.undoredo.insertOperation(unsetAtt)

        if attribute['applyOnFace']:
            patches = self.hemodel.getPatches()
            for patch in patches:
                if patch.isSelected() and not patch.isDeleted:
                    patch.setSelected(False)
                    if attribute in patch.attributes:
                        unsetAtt = UnSetAttribute(patch, attribute)
                        unsetAtt.execute()
                        self.undoredo.insertOperation(unsetAtt)

        self.undoredo.end()
        self.isChanged = True

    def addAttribute(self, _prototype, _name):
        check = self.attManager.createAttributeFromPrototype(_prototype, _name)
        if check:
            self.isChanged = True
        return check

    def saveAtribute(self, _name, _values):
        self.attManager.setAttributeValues(_name, _values)
        self.isChanged = True

    def removeAttribute(self, _name):
        self.undoredo.begin()
        delAttribute = DelAttribute(
            self.attManager, _name, self.hemodel)
        delAttribute.execute()
        self.undoredo.insertOperation(delAttribute)
        self.undoredo.end()
        self.isChanged = True

    def renameAttribute(self, _oldname, _newname):
        attributes = self.attManager.getAttributes()

        for att in attributes:
            if att['name'] == _newname:
                return False

        attribute = self.attManager.getAttributeByName(_oldname)
        attribute['name'] = _newname
        self.isChanged = True
        return True

    def setNumberOfSubdivisions(self, _number, _ratio):
        nsudv_dict = {
            "type": "Number of Subdivisions",
            "symbol": "Nsbdvs",
            "name": "Nsbdvs",
            "properties": {
                "Value": _number,
                "Ratio": _ratio,
                "Color": [0, 0, 0]
            },
            "properties_type": ["int", "float", "color"],
            "applyOnVertex": False,
            "applyOnEdge": True,
            "applyOnFace": False
        }

        self.undoredo.begin()
        segments = self.hemodel.getSegments()

        for seg in segments:
            if seg.isSelected():

                setNumber = SetNumberOfSubdivisions(seg, nsudv_dict)
                setNumber.execute()
                self.undoredo.insertOperation(setNumber)

                setAtt = SetAttribute(seg, nsudv_dict)
                setAtt.execute()
                self.undoredo.insertOperation(setAtt)

                # update mesh
                face1 = seg.edge.he1.loop.face
                face2 = seg.edge.he2.loop.face

                if face1.patch.mesh is not None:
                    self.delMesh(face1)

                if face2.patch.mesh is not None:
                    self.delMesh(face2)

        self.undoredo.end()
        self.isChanged = True

    def getAttributeSymbol(self, _attribute, _scale, _pt=None, _seg=None, _patch=None):
        return AttribSymbols.getSymbol(_attribute, _scale, _pt, _seg, _patch)


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ------------------------------- HEVIEW CLASS -------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


class HeView:

    def __init__(self, _hemodel):
        self.hemodel = _hemodel

    def getPoints(self):
        return self.hemodel.getPoints()

    def getSegments(self):
        return self.hemodel.getSegments()

    def getPatches(self):
        return self.hemodel.getPatches()

    def isEmpty(self):
        return self.hemodel.isEmpty()

    def getSelectedPoints(self):
        points = self.hemodel.points
        selected_points = []
        for pt in points:
            if pt.isSelected():
                selected_points.append(pt)

        return selected_points

    def getSelectedSegments(self):
        segments = self.hemodel.segments
        selected_segments = []
        for seg in segments:
            if seg.isSelected():
                selected_segments.append(seg)

        return selected_segments

    def getSelectedPatches(self):
        patches = self.hemodel.patches
        selected_patches = []
        for patch in patches:
            if patch.isSelected():
                selected_patches.append(patch)

        return selected_patches

    def getEntityAttributes(self, _entity):
        return _entity.attributes

    def getMeshPoints(self, _patch):
        if _patch.mesh is not None:
            return _patch.mesh.model.getPoints()

    def getMeshSegments(self, _patch):
        if _patch.mesh is not None:
            return _patch.mesh.model.getSegments()

    def getMeshPatches(self, _patch):
        if _patch.mesh is not None:
            return _patch.mesh.model.getPatches()

    def getBoundBox(self):

        if self.hemodel.isEmpty():
            return 0.0, 10.0, 0.0, 10.0

        points = self.hemodel.points
        x = points[0].getX()
        y = points[0].getY()

        xmin = x
        ymin = y
        xmax = x
        ymax = y

        for i in range(1, len(points)):
            x = points[i].getX()
            y = points[i].getY()
            xmin = min(x, xmin)
            xmax = max(x, xmax)
            ymin = min(y, ymin)
            ymax = max(y, ymax)

        for segment in self.hemodel.segments:
            xmin_c, xmax_c, ymin_c, ymax_c = segment.getBoundBox()
            xmin = min(xmin_c, xmin)
            xmax = max(xmax_c, xmax)
            ymin = min(ymin_c, ymin)
            ymax = max(ymax_c, ymax)

        return xmin, xmax, ymin, ymax

    def snapToSegment(self, _x, _y, _tol):

        if self.isEmpty():
            return False, _x, _y

        xClst = _x
        yClst = _y
        id_target = -1
        dmin = _tol

        for i in range(0, len(self.hemodel.segments)):
            xC, yC, dist = self.hemodel.segments[i].closestPoint(_x, _y)
            if dist < dmin:
                xClst = xC
                yClst = yC
                dmin = dist
                id_target = i

        if id_target < 0:
            return False, xClst, yClst

        # try to attract to a corner of the segment
        seg_pts = self.hemodel.segments[id_target].getPoints()

        dmin = _tol*2
        for pt in seg_pts:
            pt_x = pt.getX()
            pt_y = pt.getY()
            d = math.sqrt((_x-pt_x)*(_x-pt_x)+(_y-pt_y)*(_y-pt_y))

            if d < dmin:
                xClst = pt_x
                yClst = pt_y
                dmin = d

        # If found a closest point, return its coordinates
        return True, xClst, yClst

    def snapToPoint(self, _x, _y, _tol):
        if self.isEmpty():
            return False, _x, _y

        xClst = _x
        yClst = _y
        id_target = -1
        dmin = _tol

        points = self.hemodel.points
        for i in range(0, len(points)):
            xC = points[i].getX()
            yC = points[i].getY()
            if (abs(_x - xC) < _tol) and (abs(_y - yC) < _tol):
                d = math.sqrt((_x-xC)*(_x-xC)+(_y-yC)*(_y-yC))
                if d < dmin:
                    xClst = xC
                    yClst = yC
                    dmin = d
                    id_target = i

        if id_target < 0:
            return False, xClst, yClst

        # If found a closest point, return its coordinates
        return True, xClst, yClst

    def getIncidentSegmentsFromPoint(self, _point):
        incidentEdges = _point.vertex.incidentEdges()
        incidentSegments = []

        for edge in incidentEdges:
            incidentSegments.append(edge.segment)

        return incidentSegments

    def getIncidentPatchesFromPoint(self, _point):
        incidentFaces = _point.vertex.incidentFaces()
        incidentPatches = []

        for face in incidentFaces:
            if len(face.patch.segments) > 0:
                incidentPatches.append(face.patch)

        return incidentPatches

    def getAdjacentPointsFromPoint(self, _point):
        adjacentVertices = _point.vertex.adjacentVertices()
        adjacentPoints = []

        for vertex in adjacentVertices:
            adjacentPoints.append(vertex.point)

        return adjacentPoints

    def getAdjacentSegmentsFromSegment(self, _segment):
        adjacentEdges = _segment.edge.adjacentEdges()
        adjacentSegments = []

        for edge in adjacentEdges:
            adjacentSegments.append(edge.segment)

        return adjacentSegments

    def getIncidentPatchesFromSegment(self, _segment):
        incidentFaces = _segment.edge.incidentFaces()
        adjacentPatches = []

        for face in incidentFaces:
            if len(face.patch.segments) > 0:
                adjacentPatches.append(face.patch)

        return adjacentPatches

    def getIncidentPointsFromSegment(self, _segment):
        incidentVertices = _segment.edge.incidentVertices()
        adjacentPoints = []

        for vertex in incidentVertices:
            adjacentPoints.append(vertex.point)

        return adjacentPoints

    def getIncidentSegmentsFromPatch(self, _patch):
        incidentEdges = _patch.face.incidentEdges()
        adjacentSegments = []

        for edge in incidentEdges:
            adjacentSegments.append(edge.segment)

        return adjacentSegments

    def getAdjacentPatchesFromPatch(self, _patch):
        adjacentFaces = _patch.face.adjacentFaces()
        adjacentPatches = []

        for face in adjacentFaces:
            if len(face.patch.segments) > 0:
                adjacentPatches.append(face.patch)

        return adjacentPatches

    def getIncidentPointsFromPatch(self, _patch):
        incidentVertices = _patch.face.incidentVertices()
        adjacentPoints = []

        for vertex in incidentVertices:
            adjacentPoints.append(vertex.point)

        return adjacentPoints

    def getInternalPacthesFromPatch(self, _patch):
        internalFaces = _patch.face.internalFaces()
        internalPatches = []

        for face in internalFaces:
            internalPatches.append(face.patch)

        return internalPatches
