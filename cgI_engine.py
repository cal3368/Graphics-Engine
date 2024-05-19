import buffer
import math
import numpy as np
from bitarray import bitarray
import glm
import PIL

from vertex import *
from rit_window import *


class CGIengine:
    def __init__(self, myWindow, defaction):
        self.w_width = myWindow.width
        self.w_height = myWindow.height
        self.right = self.w_width
        self.left = 0
        self.bottom = 0
        self.top = self.w_height
        self.xv_min = 0
        self.yv_min = 0
        self.depth = [[1 for _ in range(self.w_height)] for _ in range(self.w_width)]
        self.image = [[[0, 0, 0] for _ in range(self.w_height)] for _ in range(self.w_width)]
        self.win = myWindow
        self.keypressed = 1
        self.default_action = defaction
        
    # go is called on every update of the window display loop
    # have your engine draw stuff in the window.
    def go(self):
        if (self.keypressed == 1):
            # default scene
            self.default_action()
        
        if (self.keypressed == 2):
            # clear the framebuffer
            self.win.clearFB (0, 0, 0)

        # push the window's framebuffer to the window
        self.win.applyFB()
        
    def keyboard (self, key) :
        if (key == '1'):
            self.keypressed = 1
            self.go()
        if (key == '2'):
            self.keypressed = 2
            self.go()

    def clearFB (self, r, g, b):
        for x in range(self.w_width):
            for y in range(self.w_height):
                self.win.set_pixel(x, y, r, g, b)

    def rasterizeLine(self, x0, y0, x1, y1, r, g, b):
        x0 = int(x0)
        x1 = int(x1)
        y0 = int(y0)
        y1 = int(y1)
        deltaX = x1 - x0
        deltaY = y1 - y0
        if deltaX != 0:
            slope = deltaY/deltaX
        else:
            slope = None
        if deltaY == 0:
            for x in range(x0, x1 + np.sign(deltaX), np.sign(deltaX)):
                self.win.set_pixel(x, y0, r, g, b)
        elif deltaX == 0:
            for y in range(y0, y1 + np.sign(deltaY), np.sign(deltaY)):
                self.win.set_pixel(x0, y, r, g, b)
        elif 0 < slope <= 1:
            incE = 2 * abs(deltaY)
            incNE = 2 * (abs(deltaY) - abs(deltaX))
            d = 2 * abs(deltaY) - abs(deltaX)
            y = y0
            for x in range(x0, x1 + np.sign(deltaX), np.sign(deltaX)):
                self.win.set_pixel(x, y, r, g, b)
                if d <= 0:
                    d += incE
                else:
                    y += np.sign(deltaY)
                    d += incNE
        elif slope > 1:
            x = x0
            incE = 2 * abs(deltaX)
            incNE = 2 * (abs(deltaX) - abs(deltaY))
            d = 2 * abs(deltaX) - abs(deltaY)
            for y in range(y0, y1 + np.sign(deltaY), np.sign(deltaY)):
                self.win.set_pixel(x, y, r, g, b)
                if d <= 0:
                    d += incE
                else:
                    d += incNE
                    x += np.sign(deltaX)
        elif -1 <= slope < 0:
            y = y0
            incE = 2 * deltaY
            incSE = 2 * (deltaY + deltaX)
            d = 2 * deltaY - deltaX
            for x in range(x0, x1 + np.sign(deltaX), np.sign(deltaX)):
                self.win.set_pixel(x, y, r, g, b)
                if d < 0:
                    d += incSE
                    y += np.sign(deltaY)
                else:
                    d += incE
        else:
            x = x0
            incE = 2 * abs(deltaX)
            incNE = 2 * (abs(deltaX) - abs(deltaY))
            d = 2 * abs(deltaX) - abs(deltaY)
            for y in range(y0, y1 + np.sign(deltaY), np.sign(deltaY)):
                self.win.set_pixel(x, y, r, g, b)
                if d <= 0:
                    d += incE
                else:
                    d += incNE
                    x += np.sign(deltaX)

    def rasterizeInterpolatedLine(self, x0, x1, y0, y1, r0, g0, b0, r1, g1, b1):
        deltaX = x1 - x0
        deltaY = y1 - y0
        deltaR = r1 - r0
        deltaG = g1 - g0
        deltaB = b1 - b0
        if deltaX != 0:
            slope = deltaY / deltaX
        else:
            slope = None
        if deltaY == 0:
            uR = deltaR / deltaX
            uG = deltaG / deltaX
            uB = deltaB / deltaX
            uR1 = deltaR / deltaX
            uG1 = deltaG / deltaX
            uB1 = deltaB / deltaX
            for x in range(x0, x1 + np.sign(deltaX), np.sign(deltaX)):
                self.win.set_pixel(x, y0, r0 + uR1, g0 + uG1, b0 + uB1)
                uR1 += uR
                uG1 += uG
                uB1 += uB
        elif deltaX == 0:
            uR = deltaR / deltaY
            uG = deltaG / deltaY
            uB = deltaB / deltaY
            uR1 = deltaR / deltaY
            uG1 = deltaG / deltaY
            uB1 = deltaB / deltaY
            for y in range(y0, y1 + np.sign(deltaY), np.sign(deltaY)):
                self.win.set_pixel(x0, y, r0 + uR1, g0 + uG1, b0 + uB1)
                uR1 += uR
                uG1 += uG
                uB1 += uB
        elif 0 < slope <= 1:
            incE = 2 * abs(deltaY)
            incNE = 2 * (abs(deltaY) - abs(deltaX))
            d = 2 * abs(deltaY) - abs(deltaX)
            y = y0
            uR = deltaR / deltaX
            uG = deltaG / deltaX
            uB = deltaB / deltaX
            uR1 = deltaR / deltaX
            uG1 = deltaG / deltaX
            uB1 = deltaB / deltaX
            for x in range(x0, x1 + np.sign(deltaX), np.sign(deltaX)):
                self.win.set_pixel(x, y, r0 + uR1, g0 + uG1, b0 + uB1)
                if d <= 0:
                    d += incE
                else:
                    y += np.sign(deltaY)
                    d += incNE
                uR1 += uR
                uG1 += uG
                uB1 += uB
        elif slope > 1:
            x = x0
            incE = 2 * abs(deltaX)
            incNE = 2 * (abs(deltaX) - abs(deltaY))
            d = 2 * abs(deltaX) - abs(deltaY)
            uR = deltaR / deltaY
            uG = deltaG / deltaY
            uB = deltaB / deltaY
            uR1 = deltaR / deltaY
            uG1 = deltaG / deltaY
            uB1 = deltaB / deltaY
            for y in range(y0, y1 + np.sign(deltaY), np.sign(deltaY)):
                self.win.set_pixel(x, y, r0 + uR1, g0 + uG1, b0 + uB1)
                if d <= 0:
                    d += incE
                else:
                    d += incNE
                    x += np.sign(deltaX)
                uR1 += uR
                uG1 += uG
                uB1 += uB
        elif -1 <= slope < 0:
            y = y0
            incE = 2 * deltaY
            incSE = 2 * (deltaY + deltaX)
            d = 2 * deltaY - deltaX
            uR = deltaR / deltaX
            uG = deltaG / deltaX
            uB = deltaB / deltaX
            uR1 = deltaR / deltaX
            uG1 = deltaG / deltaX
            uB1 = deltaB / deltaX
            for x in range(x0, x1 + np.sign(deltaX), np.sign(deltaX)):
                self.win.set_pixel(x, y, r0 + uR1, g0 + uG1, b0 + uB1)
                if d < 0:
                    d += incSE
                    y += np.sign(deltaY)
                else:
                    d += incE
                uR1 += uR
                uG1 += uG
                uB1 += uB
        else:
            x = x0
            incE = 2 * abs(deltaX)
            incNE = 2 * (abs(deltaX) - abs(deltaY))
            d = 2 * abs(deltaX) - abs(deltaY)
            uR = deltaR / deltaY
            uG = deltaG / deltaY
            uB = deltaB / deltaY
            uR1 = deltaR / deltaY
            uG1 = deltaG / deltaY
            uB1 = deltaB / deltaY
            for y in range(y0, y1 + np.sign(deltaY), np.sign(deltaY)):
                self.win.set_pixel(x, y, r0 + uR1, g0 + uG1, b0 + uB1)
                if d <= 0:
                    d += incE
                else:
                    d += incNE
                    x += np.sign(deltaX)
                uR1 += uR
                uG1 += uG
                uB1 += uB

    def rasterizeInterpolatedLines(self, vertices, colors, n):
        x0 = y0 = x1 = y1 = r0 = g0 = b0 = r1 = g1 = b1 = None
        while n > 0:
            for i in range(4):
                x0 = vertices[0]
                y0 = vertices[1]
                x1 = vertices[2]
                y1 = vertices[3]
            for i in range(6):
                r0 = colors[0]
                g0 = colors[1]
                b0 = colors[2]
                r1 = colors[3]
                g1 = colors[4]
                b1 = colors[5]
            self.rasterizeInterpolatedLine(x0, x1, y0, y1, r0, g0, b0, r1, g1, b1)
            n -= 1
            vertices = vertices[4:]
            colors = colors[6:]

    def edgeFunction(self, x, y, edge):
        return ((x - edge[0][0]) * (edge[1][1] - edge[0][1]) - (y - edge[0][1])
                * (edge[1][0] - edge[0][0]))

    def calculateArea(self, x, y, edge):
        return 1/2 * self.edgeFunction(x, y, edge)

    def calculateBarycentricCoordinates(self, x, y, edges):
        area = self.calculateArea(edges[0][0][0], edges[0][0][1], edges[1])
        if area == 0:
            area = 0.0001
        lambda0 = self.calculateArea(x, y, edges[1]) / area
        lambda1 = self.calculateArea(x, y, edges[2]) / area
        lambda2 = self.calculateArea(x, y, edges[0]) / area
        return lambda0, lambda1, lambda2

    def calculateColor(self, x, y, edges, p0, p1, p2):
        lambda0, lambda1, lambda2 = (
            self.calculateBarycentricCoordinates(x, y, edges))
        red = lambda0 * p0.r + lambda1 * p1.r + lambda2 * p2.r
        green = lambda0 * p0.g + lambda1 * p1.g + lambda2 * p2.g
        blue = lambda0 * p0.b + lambda1 * p1.b + lambda2 * p2.b
        return red, green, blue

    def rasterizeTriangle(self, p0, p1, p2):
        xMin = int(min(p0.x, p1.x, p2.x))
        xMax = math.ceil(max(p0.x, p1.x, p2.x))
        yMin = int(min(p0.y, p1.y, p2.y))
        yMax = math.ceil(max(p0.y, p1.y, p2.y))
        edges = [[(p0.x, p0.y), (p1.x, p1.y)], [(p1.x, p1.y), (p2.x, p2.y)],
                 [(p2.x, p2.y), (p0.x, p0.y)]]
        for x1 in range(xMin, xMax):
            for y in range(yMin, yMax):
                inside = True
                for edge in edges:
                    if self.edgeFunction(x1, y, edge) < 0:
                        inside = False
                        break
                if inside:
                    r, g, b = self.calculateColor(x1, y, edges, p0, p1, p2)
                    self.win.set_pixel(x1, y, r, g, b)

    def identity(self):
        return glm.identity(glm.mat4)

    def translate(self, x, y):
        return glm.translate(glm.mat3(1.0), glm.vec2(x, y))

    def scale(self, x, y):
        return glm.scale(glm.mat3(1.0), glm.vec2(x, y))

    def rotate(self, angle):
        return glm.rotate(glm.mat3(1.0), glm.radians(angle))

    def normalize(self, t, b, r, l):
        return glm.mat3x3(2/(r-l), 0, 0, 0, 2/(t-b), 0, -2*l/(r-l)-1,
                          -2*b/(t-b)-1, 1)

    def defineViewWindow(self, t, b, r, l):
        self.top = t
        self.bottom = b
        self.right = r
        self.left = l
        self.w_width = r - l
        self.w_height = t - b

    def computeOutcode(self, v, t, b, r, l):
        outcode = bitarray('0000')
        if v.y > t:
            outcode[0] = 1
        if v.y < b:
            outcode[1] = 1
        if v.x > r:
            outcode[2] = 1
        if v.x < l:
            outcode[3] = 1
        return outcode

    def clipLine(self, P0, P1, top, bottom, right, left):
        if (P1.x - P0.x) == 0:
            m = None
            B = None
        else:
            m = (P1.y - P0.y) / (P1.x - P0.x)
            B = P0.y - m * P0.x
        vertex0 = self.computeOutcode(P0, top, bottom, right, left)
        vertex1 = self.computeOutcode(P1, top, bottom, right, left)
        if (vertex0 & vertex1) != bitarray('0000'):
            return []
        if (vertex0 | vertex1) == bitarray('0000'):
            return [P0, P1]
        if vertex0 != bitarray('0000'):
            if vertex0[0] == 1:
                if not m:
                    P0.x = P0.x
                else:
                    P0.x = (top - B) / m
                P0.y = top
            if vertex0[1] == 1:
                P0.y = bottom
                if not m:
                    P0.x = P0.x
                else:
                    P0.x = (bottom - B) / m
            if vertex0[2] == 1:
                P0.x = right
                P0.y = m * right + B
            if vertex0[3] == 1:
                P0.x = left
                P0.y = m * left + B
        if vertex1 != bitarray('0000'):
            if vertex1[0] == 1:
                P1.y = top
                if not m:
                    P1.x = P1.x
                else:
                    P1.x = (top - B) / m
            if vertex1[1] == 1:
                P1.y = bottom
                if not m:
                    P1.x = P1.x
                else:
                    P1.x = (bottom - B) / m
            if vertex1[2] == 1:
                P1.x = right
                P1.y = m * right + B
            if vertex1[3] == 1:
                P1.x = left
                P1.y = m * left + B
        return [P0, P1]

    def clipPoly(self, vertices, top, bottom, right, left):
        retPoly = self.shpc(vertices, 0, top)
        if len(retPoly) > 0:
            retPoly = self.shpc(retPoly, 1, bottom)
        if len(retPoly) > 0:
            retPoly = self.shpc(retPoly, 2, right)
        if len(retPoly) > 0:
            retPoly = self.shpc(retPoly, 3, left)
        return self.polyToTriangles(retPoly)

    def shpc(self, vertices, n, side):
        vertices.append(vertices[0])
        retPoly = []
        if n == 0:
            idx = 0
            while idx < len(vertices) - 1:
                v0 = vertices[idx]
                v1 = vertices[idx+1]
                if v0.y <= side and v1.y <= side:
                    retPoly.append(v1)
                elif v0.y > side and v1.y > side:
                    pass
                elif v1.y > side:
                    newVertex = self.calculateIntersection(v0, v1, 0, side)
                    retPoly.append(newVertex)
                else:
                    newVertex = self.calculateIntersection(v0, v1, 0, side)
                    retPoly.append(newVertex)
                    retPoly.append(v1)
                idx += 1
        elif n == 1:
            idx = 0
            while idx < len(vertices) - 1:
                v0 = vertices[idx]
                v1 = vertices[idx + 1]
                if v0.y >= side and v1.y >= side:
                    retPoly.append(v1)
                elif v0.y < side and v1.y < side:
                    pass
                elif v1.y < side:
                    newVertex = self.calculateIntersection(v0, v1, 1, side)
                    retPoly.append(newVertex)
                else:
                    newVertex = self.calculateIntersection(v0, v1, 1, side)
                    retPoly.append(newVertex)
                    retPoly.append(v1)
                idx += 1
        elif n == 2:
            idx = 0
            while idx < len(vertices) - 1:
                v0 = vertices[idx]
                v1 = vertices[idx + 1]
                if v0.x <= side and v1.x <= side:
                    retPoly.append(v1)
                elif v0.x > side and v1.x > side:
                    pass
                elif v1.x > side:
                    newVertex = self.calculateIntersection(v0, v1, 2, side)
                    retPoly.append(newVertex)
                else:
                    newVertex = self.calculateIntersection(v0, v1, 2, side)
                    retPoly.append(newVertex)
                    retPoly.append(v1)
                idx += 1
        else:
            idx = 0
            while idx < len(vertices) - 1:
                v0 = vertices[idx]
                v1 = vertices[idx + 1]
                if v0.x >= side and v1.x >= side:
                    retPoly.append(v1)
                elif v0.x < side and v1.x < side:
                    pass
                elif v1.x < side:
                    newVertex = self.calculateIntersection(v0, v1, 3, side)
                    retPoly.append(newVertex)
                else:
                    newVertex = self.calculateIntersection(v0, v1, 3, side)
                    retPoly.append(newVertex)
                    retPoly.append(v1)
                idx += 1
        return retPoly

    def calculateIntersection(self, v0, v1, n, side):
        if n == 0 or n == 1:
            newY = side
            if v0.x == v1.x:
                u = self.calculateU(v0, v1, v0.x, newY)
                newR = (1-u) * v0.r + u * v1.r
                newG = (1 - u) * v0.g + u * v1.g
                newB = (1 - u) * v0.b + u * v1.b
                return Vertex(v0.x, newY, newR, newG, newB)
            m = (v1.y - v0.y) / (v1.x - v0.x)
            B = v0.y - m * v0.x
            newX = (newY - B)/m
            u = self.calculateU(v0, v1, newX, newY)
            newR = (1 - u) * v0.r + u * v1.r
            newG = (1 - u) * v0.g + u * v1.g
            newB = (1 - u) * v0.b + u * v1.b
            return Vertex(newX, newY, newR, newG, newB)
        if n == 2 or n == 3:
            newX = side
            m = (v1.y - v0.y) / (v1.x - v0.x)
            B = v0.y - m * v0.x
            newY = m * newX + B
            u = self.calculateU(v0, v1, newX, newY)
            newR = (1 - u) * v0.r + u * v1.r
            newG = (1 - u) * v0.g + u * v1.g
            newB = (1 - u) * v0.b + u * v1.b
            return Vertex(newX, newY, newR, newG, newB)

    def calculateU(self, v0, v1, newX, newY):
        num = math.sqrt((newX - v0.x)**2 + (newY - v0.y)**2)
        den = math.sqrt((v1.x - v0.x)**2 + (v1.y - v0.y)**2)
        return num/den

    def polyToTriangles(self, vertices):
        if len(vertices) == 0:
            return []
        start = vertices[0]
        retPoly = []
        idx = 1
        while idx < len(vertices)-1:
            retPoly.extend([start, vertices[idx], vertices[idx+1]])
            idx += 1
        return retPoly

    def drawTriangles2D(self, vertex_pos, colors, indices, modelT, normT):
        transformedVertices = []
        for i in range(0, len(vertex_pos), 2):
            transformedVertices.append(modelT * glm.vec3(vertex_pos[i], vertex_pos[i+1], 1.0))
        normalizedVertices = []
        for vec in transformedVertices:
            normalizedVertices.append(normT * vec)
        viewT = glm.mat3x3(self.w_width/2, 0, 0,
                           0, self.w_height/2, 0,
                           (self.right+self.left)/2, (self.top+self.bottom)/2, 1)
        finalVertices = []
        for vec in normalizedVertices:
            finalVertices.append(viewT * vec)
        for i in range(0, len(indices), 3):
            t0 = indices[i]
            t1 = indices[i+1]
            t2 = indices[i+2]
            p0 = Vertex(finalVertices[t0][0], finalVertices[t0][1],
                               colors[t0*3], colors[t0*3+1], colors[t0*3+2])
            p1 = Vertex(finalVertices[t1][0], finalVertices[t1][1],
                               colors[t1*3], colors[t1*3+1], colors[t1*3+2])
            p2 = Vertex(finalVertices[t2][0], finalVertices[t2][1],
                               colors[t2*3], colors[t2*3+1], colors[t2*3+2])
            newVertices = self.clipPoly([p0, p1, p2], self.top, self.bottom, self.right, self.left)
            for i2 in range(0, len(newVertices), 3):
                self.rasterizeTriangle(newVertices[i2], newVertices[i2+1], newVertices[i2+2])

    def drawLines(self, vertex_pos, colors, indices, modelT, normT):
        transformedVertices = []
        for i in range(0, len(vertex_pos), 2):
            transformedVertices.append(
                modelT * glm.vec3(vertex_pos[i], vertex_pos[i + 1], 1.0))
        normalizedVertices = []
        for vec in transformedVertices:
            normalizedVertices.append(normT * vec)
        viewT = glm.mat3x3(self.w_width / 2, 0, 0, 0, self.w_height / 2, 0,
                           (self.right + self.left) / 2,
                           (self.top + self.bottom) / 2, 1)
        finalVertices = []
        for vec in normalizedVertices:
            finalVertices.append(viewT * vec)
        for i in range(0, len(indices), 2):
            t0 = indices[i]
            t1 = indices[i+1]
            p0 = Vertex(finalVertices[t0][0], finalVertices[t0][1],
                        colors[t0*2], colors[t0*2+1], colors[t0*2+2])
            p1 = Vertex(finalVertices[t1][0], finalVertices[t1][1],
                        colors[t1*2], colors[t1*2+1], colors[t1*2+2])
            newVertices = self.clipLine(p0, p1, self.top, self.bottom, self.right, self.left)
            for i2 in range(0, len(newVertices), 2):
                self.rasterizeInterpolatedLine(newVertices[i2].x, newVertices[i2+1].x,
                                               newVertices[i2].y, newVertices[i2+1].y,
                                               newVertices[i2].r, newVertices[i2].g,
                                               newVertices[i2].b, newVertices[i2+1].r,
                                               newVertices[i2+1].g, newVertices[i2+1].b)

    def identity3D(self):
        return glm.identity(glm.mat4)

    def scale3D(self, x, y, z):
        return glm.scale(glm.vec3(x, y, z))

    def translate3D(self, x, y, z):
        return glm.translate(glm.vec3(x, y, z))

    def rotateX(self, angle):
        return glm.rotate(glm.radians(angle), glm.vec3(1, 0, 0))

    def rotateY(self, angle):
        return glm.rotate(glm.radians(angle), glm.vec3(0, 1, 0))

    def rotateZ(self, angle):
        return glm.rotate(glm.radians(angle), glm.vec3(0, 0, 1))

    def ortho3D(self, l, r, b, t, n, f):
        return glm.orthoRH_NO(l, r, b, t, n, f)

    def lookAt(self, eye, lookat, up):
        return glm.lookAtRH(eye, lookat, up)

    def frustum3D(self, l, r, b, t, n, f):
        return glm.frustumRH_NO(l, r, b, t, n, f)

    def getDepth(self, x, y, edges, z0, z1, z2):
        l0, l1, l2 = self.calculateBarycentricCoordinates(x, y, edges)
        return l0*z0 + l1*z1 + l2*z2

    def rasterizeTriangleNew(self, p0, p1, p2, f_shader, uniforms):
        xMin = int(min(p0.x, p1.x, p2.x))
        xMax = math.ceil(max(p0.x, p1.x, p2.x))
        yMin = int(min(p0.y, p1.y, p2.y))
        yMax = math.ceil(max(p0.y, p1.y, p2.y))
        edges = [[(p0.x, p0.y), (p1.x, p1.y)], [(p1.x, p1.y), (p2.x, p2.y)],
                 [(p2.x, p2.y), (p0.x, p0.y)]]
        for x1 in range(xMin, xMax):
            for y in range(yMin, yMax):
                inside = True
                for edge in edges:
                    if self.edgeFunction(x1, y, edge) < 0:
                        inside = False
                        break
                if inside:
                    d = self.getDepth(x1, y, edges, p0.z, p1.z, p2.z)
                    if d < self.depth[x1][y]:
                        self.depth[x1][y] = d
                        l0, l1, l2 = self.calculateBarycentricCoordinates(x1, y,
                                                                          edges)
                        r, g, b = f_shader(p0, p1, p2, l0, l1, l2, uniforms)
                        self.win.set_pixel(x1, y, r, g, b)

    def createBuffer(self, name, n_per_vertex):
        return buffer.Buffer(name, n_per_vertex)

    def fillBuffer(self, id, data):
        id.data = data

    def drawTriangles(self, vertex_pos_buffer, vertex_index_buffer,
                      vertex_shader, fragment_shader, uniforms):
        posN = vertex_pos_buffer.n_per_vertex
        posData = vertex_pos_buffer.data
        indicesData = vertex_index_buffer.data
        normalsData = uniforms.get("normals").data
        uvsN = uniforms.get("uvs").n_per_vertex
        uvsData = uniforms.get("uvs").data
        verticesOriginal = []
        uvs = []
        normals = []
        for i in range(0, len(posData), posN):
            v = Vertex(posData[i], posData[i + 1], posData[i + 2])
            norm = glm.vec3(normalsData[i], normalsData[i + 1], normalsData[i + 2])
            normals.append(norm)
            verticesOriginal.append(v)
        for i in range(0, len(uvsData), uvsN):
            uvs.append(glm.vec2(uvsData[i], uvsData[i + 1]))
        verticesTransformed = []
        for i in range(len(verticesOriginal)):
            verticesOriginal[i].attachValue("normal", normals[i])
            verticesOriginal[i].attachValue("uv", uvs[i])
            verticesTransformed.append(vertex_shader(verticesOriginal[i], uniforms))
        viewPortT = glm.mat3x3(self.w_width / 2, 0, 0,
                               0, self.w_height / 2, 0,
                               (self.right + self.left) / 2,
                               (self.top + self.bottom) / 2, 1)
        for v in verticesTransformed:
            v.x = v.x / v.getValue("W")
            v.y = v.y / v.getValue("W")
            v.z = v.z / v.getValue("W")
            newVec = glm.vec3(1)
            newVec.x = v.x
            newVec.y = v.y
            final = viewPortT * newVec
            v.x = final.x
            v.y = final.y
        for i in range(0, len(indicesData), 3):
            t0 = indicesData[i]
            t1 = indicesData[i + 1]
            t2 = indicesData[i + 2]
            p0 = verticesTransformed[t0]
            p1 = verticesTransformed[t1]
            p2 = verticesTransformed[t2]
            self.rasterizeTriangleNew(p0, p1, p2, fragment_shader, uniforms)

