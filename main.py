from rit_window import *
from cgI_engine import *
from vertex import *
# from clipper import *
from shapes import *
import numpy as np
from PIL import Image


baseColor = glm.vec3(153/255, 170/255, 186/255)
legColor = glm.vec3(80/255, 86/255, 96/255)


def vertexShader(V, uniforms):
    modelT = uniforms.get("modelT")
    viewT = uniforms.get("viewT")
    projectionT = uniforms.get("projectionT")
    camCoords = viewT * modelT * glm.vec4(V.x, V.y, V.z, 1.0)
    V.attachValue("Coordinates", camCoords)
    proj = projectionT * camCoords
    V.values["W"] = proj.w
    V.x = proj.x
    V.y = proj.y
    V.z = proj.z
    return V


def imageShader(p0, p1, p2, alpha, beta, gamma, uniforms):
    im = uniforms.get("image")
    u0, v0 = p0.values.get("uv")[0], p0.values.get("uv")[1]
    u1, v1 = p1.values.get("uv")[0], p1.values.get("uv")[1]
    u2, v2 = p2.values.get("uv")[0], p2.values.get("uv")[1]
    newU = (alpha * u0 + beta * u1 + gamma * u2) * (im.width - 1)
    newV = (alpha * v0 + beta * v1 + gamma * v2) * (im.height - 1)
    return (im.getpixel((newU, newV))[0]/255,
               im.getpixel((newU, newV))[1]/255,
               im.getpixel((newU, newV))[2]/255)


def phongShader(p0, p1, p2, alpha, beta, gamma, uniforms):
    nX0, nY0, nZ0 = p0.values.get("normal")[0], p0.values.get("normal")[1], \
    p0.values.get("normal")[2]
    nX1, nY1, nZ1 = p1.values.get("normal")[0], p1.values.get("normal")[1], \
    p1.values.get("normal")[2]
    nX2, nY2, nZ2 = p2.values.get("normal")[0], p2.values.get("normal")[1], \
    p2.values.get("normal")[2]
    oldX0, oldY0, oldZ0 = p0.values.get("Coordinates").x, p0.values.get(
        "Coordinates").y, p0.values.get("Coordinates").z
    oldX1, oldY1, oldZ1 = p1.values.get("Coordinates").x, p1.values.get(
        "Coordinates").y, p1.values.get("Coordinates").z
    oldX2, oldY2, oldZ2 = p2.values.get("Coordinates").x, p2.values.get(
        "Coordinates").y, p2.values.get("Coordinates").z
    newNX = alpha * nX0 + beta * nX1 + gamma * nX2
    newNY = alpha * nY0 + beta * nY1 + gamma * nY2
    newNZ = alpha * nZ0 + beta * nZ1 + gamma * nZ2
    newX = alpha * oldX0 + beta * oldX1 + gamma * oldX2
    newY = alpha * oldY0 + beta * oldY1 + gamma * oldY2
    newZ = alpha * oldZ0 + beta * oldZ1 + gamma * oldZ2
    newNorm = glm.vec3(newNX, newNY, newNZ)
    vec = glm.vec3(newX, newY, newZ)
    kA, kD, kS = uniforms.get("k")[0], uniforms.get("k")[1], uniforms.get("k")[2]
    oR, oG, oB = uniforms.get("ocolor")[0], uniforms.get("ocolor")[1], uniforms.get("ocolor")[2]
    sR, sG, sB = uniforms.get("scolor")[0], uniforms.get("scolor")[1], uniforms.get("scolor")[2]
    lR, lG, lB = uniforms.get("lightcolor")[0], uniforms.get("lightcolor")[1], uniforms.get("lightcolor")[2]
    lX, lY, lZ = uniforms.get("lightpos").x, uniforms.get("lightpos").y, uniforms.get("lightpos").z
    ambR, ambG, ambB = uniforms.get("amb_color")[0], uniforms.get("amb_color")[1], uniforms.get("amb_color")[2]
    AR, AG, AB = [oR * ambR, oG * ambG, oB * ambB]
    lVec = glm.vec3(lX - vec.x, lY - vec.y, lZ - vec.z)
    M = glm.mat4x4(0.5, 0, 0, 0,
                   0, 1, 0, 0,
                   0, 0, 1, 0,
                   0, 0, 0, 1)
    normal = M * newNorm
    nDotL = glm.dot(glm.normalize(normal), glm.normalize(lVec))
    if nDotL < 0:
        nDotL = 0
    DR, DG, DB = [oR * lR * nDotL, oG * lG * nDotL, oB * lB * nDotL]
    vVec = glm.vec3(0 - vec.x, 0 - vec.y, 0 - vec.z)
    rVec = glm.reflect(glm.normalize(lVec * -1), glm.normalize(normal))
    spec = glm.dot(glm.normalize(rVec), glm.normalize(vVec))
    if spec < 0:
        spec = 0
    spec = spec ** uniforms.get("exponent")
    SR, SG, SB = [lR * sR * spec, lG * sG * spec, lB * sB * spec]
    return [kA * AR + kD * DR + kS * SR, kA * AG + kD * DG + kS * SG,
            kA * AB + kD * DB + kS * SB]


def drawHead(viewT, projectionT):
    uniform = {}
    uniform["viewT"] = viewT
    uniform["projectionT"] = projectionT
    uniform["ocolor"] = baseColor
    uniform["scolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["k"] = glm.vec3(0.6, 0.4, 0)
    uniform["exponent"] = 1
    uniform["lightpos"] = glm.vec3(0.0, 5.0, -2.0)
    uniform["lightcolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["amb_color"] = glm.vec3(1.0, 1.0, 1.0)

    uniform["normals"] = normalSphereBuffer
    uniform["uvs"] = uvSphereBuffer

    modelT = myEngine.translate3D(0, 3.5, -10) * myEngine.scale3D(0.05, 0.05, 0.05) * myEngine.rotateY(90)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(sphereBuffer, indexSphereBuffer, vertexShader, phongShader,
                           uniform)

    uniform["normals"] = normalConeBuffer
    uniform["uvs"] = uvConeBuffer

    modelT = myEngine.translate3D(0, 3.4, -10.0) * myEngine.scale3D(0.05, 0.25, 0.05) * myEngine.rotateY(90)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(coneBuffer, coneIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalSphereBuffer
    uniform["uvs"] = uvSphereBuffer

    modelT = myEngine.translate3D(0, 3.25, -10.0) * myEngine.scale3D(0.12, 0.05,
                                                                   0.12) * myEngine.rotateY(90)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(sphereBuffer, indexSphereBuffer, vertexShader,
                           phongShader,
                           uniform)

    modelT = myEngine.translate3D(0, 2.92, -10.0) * myEngine.scale3D(0.60, 0.65,
                                                                    0.60) * myEngine.rotateY(90)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(sphereBuffer, indexSphereBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalCylinderBuffer
    uniform["uvs"] = uvCylinderBuffer

    modelT = myEngine.translate3D(0, 2.57, -10.0) * myEngine.scale3D(0.60, 0.7,
                                                                    0.60) * myEngine.rotateY(90)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalCubeBuffer
    uniform["uvs"] = uvCubeBuffer

    modelT = myEngine.translate3D(-0.1, 2.55, -8.0) * myEngine.scale3D(0.4, 0.2, 0) * myEngine.rotateX(-10) * myEngine.rotateY(-20) * myEngine.rotateZ(1)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cubeBuffer, cubeIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["image"] = Image.open("benderEyes.png")
    modelT = myEngine.translate3D(-0.19, 2.56, -7.99) * myEngine.rotateZ(3) * myEngine.rotateY(-3.4) * myEngine.rotateZ(1) * myEngine.scale3D(0.4, 0.2, 0)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cubeBuffer, cubeIndexBuffer, vertexShader,
                           imageShader,
                           uniform)

    uniform["image"] = Image.open("benderMouth.png")
    uniform["normals"] = normalCylinderBuffer
    uniform["uvs"] = uvCylinderBuffer

    modelT = myEngine.translate3D(0, 2.4, -9.99) * myEngine.scale3D(0.60, 0.4,
                                                                     0.60) * myEngine.rotateY(
        58)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           imageShader,
                           uniform)


def drawBody(viewT, projectionT):
    uniform = {}
    uniform["viewT"] = viewT
    uniform["projectionT"] = projectionT
    uniform["ocolor"] = baseColor
    uniform["scolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["k"] = glm.vec3(0.6, 0.4, 0)
    uniform["exponent"] = 1
    uniform["lightpos"] = glm.vec3(0.0, 5.0, -2.0)
    uniform["lightcolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["amb_color"] = glm.vec3(1.0, 1.0, 1.0)

    uniform["normals"] = normalConeBuffer
    uniform["uvs"] = uvConeBuffer

    modelT = myEngine.translate3D(0, 2.2, -10.0) * myEngine.scale3D(1.3,0.4,1.3)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(coneBuffer, coneIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalCylinderBuffer
    uniform["uvs"] = uvCylinderBuffer

    modelT = myEngine.translate3D(0, 0.90, -10.0) * myEngine.scale3D(1.3, 2.2,
                                                                     1.3) * myEngine.rotateY(65)
    uniform["modelT"] = modelT
    uniform["image"] = Image.open("benderBody.png")
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           imageShader,
                           uniform)

    modelT = myEngine.translate3D(0, 2, -10.0) * myEngine.scale3D(1.3, 0.015,
                                                                     1.3)
    uniform["modelT"] = modelT
    uniform["ocolor"] = glm.vec3(0, 0, 0)
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    modelT = myEngine.translate3D(0, 2.2, -10) * myEngine.scale3D(0.6, 0.01,
                                                                  0.6)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalSphereBuffer
    uniform["uvs"] = uvSphereBuffer
    uniform["ocolor"] = baseColor

    modelT = myEngine.translate3D(0.55, 1.76, -9.6) * myEngine.scale3D(0.35, 0.45, 0.45)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(sphereBuffer, indexSphereBuffer, vertexShader,
                           phongShader,
                           uniform)

    modelT = myEngine.translate3D(-0.6, 1.76, -10.4) * myEngine.scale3D(0.35, 0.45,
                                                                      0.45)
    uniform["modelT"] = modelT
    myEngine.drawTriangles(sphereBuffer, indexSphereBuffer, vertexShader,
                           phongShader,
                           uniform)


def drawLeftLeg(viewT, projectionT):
    uniform = {}
    uniform["viewT"] = viewT
    uniform["projectionT"] = projectionT
    uniform["ocolor"] = legColor
    uniform["scolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["k"] = glm.vec3(0.6, 0.4, 0)
    uniform["exponent"] = 1
    uniform["lightpos"] = glm.vec3(0.0, 5.0, -2.0)
    uniform["lightcolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["amb_color"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["normals"] = normalCylinderBuffer
    uniform["uvs"] = uvCylinderBuffer

    modelT = (myEngine.rotateZ(8) * myEngine.translate3D(.2, -0.5, -9.8) *
              myEngine.scale3D(0.35, .5, 0.35)
              * myEngine.rotateZ(-3))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(8) * myEngine.translate3D(.2, -0.76, -9.8) *
              myEngine.scale3D(0.35, 0.03, 0.35)
              * myEngine.rotateZ(-3))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.rotateZ(8) * myEngine.translate3D(.2, -0.98, -9.8) *
              myEngine.scale3D(0.35, 0.42, 0.35)
              * myEngine.rotateZ(2))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(8) * myEngine.translate3D(.2, -1.2, -9.8) *
              myEngine.scale3D(0.35, 0.03, 0.35)
              * myEngine.rotateZ(2))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.rotateZ(8) * myEngine.translate3D(.2, -1.42, -9.8) *
              myEngine.scale3D(0.35, 0.42, 0.35)
              * myEngine.rotateZ(-2))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(8) * myEngine.translate3D(.2, -1.64, -9.8) *
              myEngine.scale3D(0.35, 0.04, 0.35)
              * myEngine.rotateZ(-2))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.rotateZ(5) * myEngine.translate3D(.279, -1.88, -9.8) *
              myEngine.scale3D(0.35, 0.5, 0.35) * myEngine.rotateZ(-2))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(5) * myEngine.translate3D(.279, -2.14, -9.8) *
              myEngine.scale3D(0.35, 0.04, 0.35)
              * myEngine.rotateZ(-2))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.rotateZ(5) * myEngine.translate3D(.27, -2.33, -9.8) *
              myEngine.scale3D(0.35, 0.35, 0.35) * myEngine.rotateZ(-2))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(5) * myEngine.translate3D(.27, -2.51, -9.8) *
              myEngine.scale3D(0.35, 0.02, 0.35)
              * myEngine.rotateZ(-2))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalSemiSphereBuffer
    uniform["uvs"] = uvSemiSphereBuffer
    uniform["ocolor"] = baseColor

    modelT = (myEngine.rotateZ(5) * myEngine.translate3D(.27, -2.8, -9.8) *
              myEngine.scale3D(0.9, 0.5, 0.9)
              * myEngine.rotateZ(-2) * myEngine.rotateX(-20))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(semiSphereBuffer, indexSemiSphereBuffer, vertexShader,
                           phongShader,
                           uniform)


def drawRightLeg(viewT, projectionT):
    uniform = {}
    uniform["viewT"] = viewT
    uniform["projectionT"] = projectionT
    uniform["ocolor"] = legColor
    uniform["scolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["k"] = glm.vec3(0.6, 0.4, 0)
    uniform["exponent"] = 1
    uniform["lightpos"] = glm.vec3(0.0, 5.0, -2.0)
    uniform["lightcolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["amb_color"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["normals"] = normalCylinderBuffer
    uniform["uvs"] = uvCylinderBuffer

    modelT = (myEngine.rotateZ(-15) * myEngine.translate3D(-0.2, -0.62, -10.2) *
              myEngine.scale3D(0.35, .5, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(-15) * myEngine.translate3D(-0.2, -0.89, -10.2) *
              myEngine.scale3D(0.35, 0.02, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.rotateZ(-14) * myEngine.translate3D(-0.22, -1.14, -10.2) *
              myEngine.scale3D(0.35, 0.5, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(-14) * myEngine.translate3D(-.22, -1.4, -10.2) *
              myEngine.scale3D(0.35, 0.03, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.rotateZ(-8) * myEngine.translate3D(-.36, -1.6, -10.2) *
              myEngine.scale3D(0.35, 0.42, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(-8) * myEngine.translate3D(-.36, -1.79, -10.2) *
              myEngine.scale3D(0.35, 0.04, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.rotateZ(-2) * myEngine.translate3D(-.55, -1.99, -10.2) *
              myEngine.scale3D(0.35, 0.47, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(-2) * myEngine.translate3D(-.55, -2.24, -10.2) *
              myEngine.scale3D(0.35, 0.04, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.rotateZ(-2) * myEngine.translate3D(-.55, -2.4, -10.2) *
              myEngine.scale3D(0.35, 0.32, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.rotateZ(-2) * myEngine.translate3D(-.55, -2.58, -10.2) *
              myEngine.scale3D(0.35, 0.02, 0.35))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalSemiSphereBuffer
    uniform["uvs"] = uvSemiSphereBuffer
    uniform["ocolor"] = baseColor

    modelT = (myEngine.rotateZ(-2) * myEngine.translate3D(-.55, -2.9, -10.2) *
              myEngine.scale3D(0.9, 0.5, 0.9) * myEngine.rotateX(-20))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(semiSphereBuffer, indexSemiSphereBuffer,
                           vertexShader,
                           phongShader,
                           uniform)


def drawLeftArm(viewT, projectionT):
    uniform = {}
    uniform["viewT"] = viewT
    uniform["projectionT"] = projectionT
    uniform["ocolor"] = legColor
    uniform["scolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["k"] = glm.vec3(0.6, 0.4, 0)
    uniform["exponent"] = 1
    uniform["lightpos"] = glm.vec3(0.0, 5.0, -2.0)
    uniform["lightcolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["amb_color"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["normals"] = normalCylinderBuffer
    uniform["uvs"] = uvCylinderBuffer

    modelT = (myEngine.translate3D(0.67, 1.68, -9.5) *
              myEngine.rotateZ(55) * myEngine.rotateX(-25) * myEngine.scale3D(0.32, .3, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(0.78, 1.6, -9.4) *
              myEngine.rotateZ(55) * myEngine.rotateX(-25) *
              myEngine.scale3D(0.32, 0.02, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.translate3D(0.93, 1.48, -9.3) *
              myEngine.rotateZ(50) * myEngine.rotateX(-25) *
              myEngine.scale3D(0.32, 0.4, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(1.07, 1.36, -9.2) *
              myEngine.rotateZ(50) * myEngine.rotateX(-25) *
              myEngine.scale3D(0.32, 0.02, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.translate3D(1.21, 1.21, -9.1) *
              myEngine.rotateZ(40) * myEngine.rotateX(-25) *
              myEngine.scale3D(0.32, 0.42, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(1.34, 1.06, -9) *
              myEngine.rotateZ(40) * myEngine.rotateX(-25) *
              myEngine.scale3D(0.32, 0.02, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.translate3D(1.43, 0.87, -8.9) *
              myEngine.rotateZ(25) * myEngine.rotateX(-25) *
              myEngine.scale3D(0.32, 0.47, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(1.51, .66, -8.8) *
              myEngine.rotateZ(20) * myEngine.rotateX(-25) *
              myEngine.scale3D(0.32, 0.02, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.translate3D(1.52, .47, -8.7) * myEngine.rotateX(-25) *
              myEngine.scale3D(0.32, 0.45, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(1.52, .26, -8.6) * myEngine.rotateX(-25) *
              myEngine.scale3D(0.32, 0.01, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalConeBuffer
    uniform["uvs"] = uvConeBuffer
    uniform["ocolor"] = baseColor

    modelT = (myEngine.translate3D(1.53, .38, -8.65) * myEngine.rotateZ(-5)
              * myEngine.rotateX(-30) *
              myEngine.scale3D(0.4, 0.6, 0.4))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(coneBuffer, coneIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalCylinderBuffer
    uniform["uvs"] = uvCylinderBuffer

    modelT = (myEngine.translate3D(1.53, 0, -8.8) *
              myEngine.scale3D(0.1, 0.4, 0.1))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    modelT = (myEngine.translate3D(1.62, 0, -8.79) * myEngine.rotateZ(5) *
              myEngine.rotateX(-10) * myEngine.scale3D(0.1, 0.4, 0.1))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    modelT = (myEngine.translate3D(1.72, 0, -8.8) * myEngine.rotateZ(15) *
              myEngine.rotateX(5) * myEngine.scale3D(0.1, 0.4, 0.1))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)


def drawRightArm(viewT, projectionT):
    uniform = {}
    uniform["viewT"] = viewT
    uniform["projectionT"] = projectionT
    uniform["ocolor"] = legColor
    uniform["scolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["k"] = glm.vec3(0.6, 0.4, 0)
    uniform["exponent"] = 1
    uniform["lightpos"] = glm.vec3(0.0, 5.0, -2.0)
    uniform["lightcolor"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["amb_color"] = glm.vec3(1.0, 1.0, 1.0)
    uniform["normals"] = normalCylinderBuffer
    uniform["uvs"] = uvCylinderBuffer

    modelT = (myEngine.translate3D(-.8, 1.65, -10.4) *
              myEngine.rotateZ(-55) * myEngine.rotateX(25) *
              myEngine.scale3D(0.32, .3, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(-0.93, 1.55, -10.5) *
              myEngine.rotateZ(-55) * myEngine.rotateX(25) *
              myEngine.scale3D(0.32, 0.02, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.translate3D(-1.11, 1.4, -10.6) *
              myEngine.rotateZ(-50) * myEngine.rotateX(25) *
              myEngine.scale3D(0.32, 0.5, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(-1.29, 1.24, -10.7) *
              myEngine.rotateZ(-50) * myEngine.rotateX(25) *
              myEngine.scale3D(0.32, 0.02, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.translate3D(-1.42, 1.02, -10.8) *
              myEngine.rotateZ(-30) * myEngine.rotateX(20) *
              myEngine.scale3D(0.32, 0.52, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(-1.55, .8, -10.9) *
              myEngine.rotateZ(-30) * myEngine.rotateX(20) *
              myEngine.scale3D(0.32, 0.02, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.translate3D(-1.62, 0.55, -11) *
              myEngine.rotateZ(-15) * myEngine.rotateX(20) *
              myEngine.scale3D(0.32, 0.52, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(-1.69, .3, -11.1) *
              myEngine.rotateZ(-15) * myEngine.rotateX(20) *
              myEngine.scale3D(0.32, 0.02, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = legColor
    modelT = (myEngine.translate3D(-1.7, 0.05, -11.2) *
              myEngine.rotateZ(2) * myEngine.rotateX(15) *
              myEngine.scale3D(0.32, 0.45, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["ocolor"] = glm.vec3(0, 0, 0)
    modelT = (myEngine.translate3D(-1.7, -.19, -11.3) *
              myEngine.rotateZ(2) * myEngine.rotateX(15) *
              myEngine.scale3D(0.32, 0.02, 0.32))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalConeBuffer
    uniform["uvs"] = uvConeBuffer
    uniform["ocolor"] = baseColor

    modelT = (myEngine.translate3D(-1.73, -.12, -11.4) * myEngine.rotateZ(2)
              * myEngine.rotateX(15) *
              myEngine.scale3D(0.4, 0.6, 0.4))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(coneBuffer, coneIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    uniform["normals"] = normalCylinderBuffer
    uniform["uvs"] = uvCylinderBuffer

    modelT = (myEngine.translate3D(-1.85, -.6, -11.6) *
              myEngine.rotateX(-10) *
              myEngine.scale3D(0.1, 0.35, 0.1))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    modelT = (myEngine.translate3D(-1.7, -.6, -11.55) * myEngine.rotateZ(-5) *
              myEngine.rotateX(10) * myEngine.scale3D(0.1, 0.35, 0.1))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)

    modelT = (myEngine.translate3D(-1.62, -.6, -11.65) * myEngine.rotateZ(15) *
              myEngine.rotateX(-5) * myEngine.scale3D(0.1, 0.4, 0.1))
    uniform["modelT"] = modelT
    myEngine.drawTriangles(cylinderBuffer, cylinderIndexBuffer, vertexShader,
                           phongShader,
                           uniform)


def default_action():
    # clear the FB
    myEngine.win.clearFB(45/255, 45/255, 45/255)

    viewT = myEngine.lookAt(glm.vec3(0.0, 2.0, -4.0), glm.vec3(0.0, 0.0, -10.0),
                            glm.vec3(0.0, 1.0, 0.0))
    projectionT = myEngine.frustum3D(-6, 6, -6, 6, 8.0, 15.0)

    drawHead(viewT, projectionT)
    drawBody(viewT, projectionT)
    drawLeftLeg(viewT, projectionT)
    drawRightLeg(viewT, projectionT)
    drawLeftArm(viewT, projectionT)
    drawRightArm(viewT, projectionT)


window = RitWindow(800, 800)
myEngine = CGIengine(window, default_action)

sphereBuffer = myEngine.createBuffer("vertices", 3)
myEngine.fillBuffer(sphereBuffer, sphere)
indexSphereBuffer = myEngine.createBuffer("indices", 1)
myEngine.fillBuffer(indexSphereBuffer, sphere_idx)
uvSphereBuffer = myEngine.createBuffer("uvs", 2)
myEngine.fillBuffer(uvSphereBuffer, sphere_uv)
normalSphereBuffer = myEngine.createBuffer("normals", 3)
myEngine.fillBuffer(normalSphereBuffer, sphere_normals)

cylinderBuffer = myEngine.createBuffer("vertices", 3)
myEngine.fillBuffer(cylinderBuffer, cylinder)
cylinderIndexBuffer = myEngine.createBuffer("indices", 1)
myEngine.fillBuffer(cylinderIndexBuffer, cylinder_idx)
uvCylinderBuffer = myEngine.createBuffer("uvs", 2)
myEngine.fillBuffer(uvCylinderBuffer, cylinder_uv)
normalCylinderBuffer = myEngine.createBuffer("normals", 3)
myEngine.fillBuffer(normalCylinderBuffer, cylinder_normals)

coneBuffer = myEngine.createBuffer("vertices", 3)
myEngine.fillBuffer(coneBuffer, cone)
coneIndexBuffer = myEngine.createBuffer("indices", 1)
myEngine.fillBuffer(coneIndexBuffer, cone_idx)
uvConeBuffer = myEngine.createBuffer("uvs", 2)
myEngine.fillBuffer(uvConeBuffer, cone_uv)
normalConeBuffer = myEngine.createBuffer("normals", 3)
myEngine.fillBuffer(normalConeBuffer, cone_normals)

cubeBuffer = myEngine.createBuffer("vertices", 3)
myEngine.fillBuffer(cubeBuffer, cube)
cubeIndexBuffer = myEngine.createBuffer("indices", 1)
myEngine.fillBuffer(cubeIndexBuffer, cube_idx)
uvCubeBuffer = myEngine.createBuffer("uvs", 2)
myEngine.fillBuffer(uvCubeBuffer, cube_uv)
normalCubeBuffer = myEngine.createBuffer("normals", 3)
myEngine.fillBuffer(normalCubeBuffer, cube_normals)

semiSphereBuffer = myEngine.createBuffer("vertices", 3)
myEngine.fillBuffer(semiSphereBuffer, sphere[:len(sphere)//2])
indexSemiSphereBuffer = myEngine.createBuffer("indices", 1)
myEngine.fillBuffer(indexSemiSphereBuffer, sphere_idx[:len(sphere_idx)//2])
uvSemiSphereBuffer = myEngine.createBuffer("uvs", 2)
myEngine.fillBuffer(uvSemiSphereBuffer, sphere_uv[:len(sphere_uv)//2])
normalSemiSphereBuffer = myEngine.createBuffer("normals", 3)
myEngine.fillBuffer(normalSemiSphereBuffer, sphere_normals[:len(sphere_normals)//2])


def main():
    window.run(myEngine)


if __name__ == "__main__":
    main()
