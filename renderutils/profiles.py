import numpy as np
from panda3d.core import GeomVertexFormat
from panda3d.core import GeomVertexData, GeomVertexWriter, Geom
from panda3d.core import GeomNode, GeomTriangles


class I_Profile():
    """
    This class is used to render the profile geometry of an I Beam profile.
    """

    def __init__(self, beam, s, base):
        """
        Initiates the I_Profile class.
        """
        self.beam = beam
        self.s = s
        self.base = base

    def create_profile(self):
        """
        Creates the necessary rendering geometry to render the 3d rendition of
        an IBeam profile.
        """
        beam = self.beam
        s = self.s
        prof = beam.profile

        T = beam.elements[0].trans_matrix[:3, :3]

        fmt = GeomVertexFormat.getV3n3c4()
        vertexData = GeomVertexData('something', fmt, Geom.UHStatic)

        vertexData.setNumRows(24)

        vertices = GeomVertexWriter(vertexData, 'vertex')
        normals = GeomVertexWriter(vertexData, 'normal')
        colors = GeomVertexWriter(vertexData, 'color')

        b_0 = np.array([beam.x[0], beam.y[0], beam.z[0]]) / s
        b_l = np.array([beam.x[-1], beam.y[-1], beam.z[-1]]) / s

        # 0
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    -prof.b / 2,
                                                    -prof.h / 2])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 1
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    prof.b / 2,
                                                    -prof.h / 2])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 2
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    prof.b / 2,
                                                    -prof.h / 2 +
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 3
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    prof.t_w/2,
                                                    -prof.h / 2 +
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
        # 4
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    prof.t_w/2,
                                                    prof.h / 2 -
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
        # 5
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    prof.b / 2,
                                                    prof.h / 2 -
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 6
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    prof.b / 2,
                                                    prof.h / 2])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 7
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    -prof.b / 2,
                                                    prof.h / 2])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 8
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    -prof.b / 2,
                                                    prof.h / 2 -
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 9
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    - prof.t_w/2,
                                                    prof.h / 2 -
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
        # 10
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    - prof.t_w/2,
                                                    - prof.h/2 +
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
        # 11
        vertices.addData3f(*(b_0 + np.dot(T.T,
                                          np.array([0,
                                                    -prof.b / 2,
                                                    - prof.h/2 +
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 12
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    -prof.b / 2,
                                                    -prof.h / 2])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 13
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    prof.b / 2,
                                                    -prof.h / 2])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 14
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    prof.b / 2,
                                                    -prof.h / 2 +
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 15
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    prof.t_w/2,
                                                    -prof.h / 2 +
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
        # 16
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    prof.t_w/2,
                                                    prof.h / 2 -
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
        # 17
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    prof.b / 2,
                                                    prof.h / 2 -
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 18
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    prof.b / 2,
                                                    prof.h / 2])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 19
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    -prof.b / 2,
                                                    prof.h / 2])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 20
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    -prof.b / 2,
                                                    prof.h / 2 -
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
        # 21
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    - prof.t_w/2,
                                                    prof.h / 2 -
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
        # 22
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    - prof.t_w/2,
                                                    - prof.h/2 +
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
        # 23
        vertices.addData3f(*(b_l + np.dot(T.T,
                                          np.array([0,
                                                    -prof.b / 2,
                                                    - prof.h/2 +
                                                    prof.t_f])
                                          / s)))
        normals.addData3f(0, 0, 1)
        colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)

        primitive = GeomTriangles(Geom.UHStatic)

        # front
        primitive.addVertices(0, 1, 2)
        primitive.addVertices(2, 1, 0)
        primitive.addVertices(0, 2, 11)
        primitive.addVertices(11, 2, 0)
        primitive.addVertices(10, 3, 4)
        primitive.addVertices(4, 3, 10)
        primitive.addVertices(10, 4, 9)
        primitive.addVertices(9, 4, 10)
        primitive.addVertices(6, 7, 5)
        primitive.addVertices(5, 7, 6)
        primitive.addVertices(7, 8, 5)
        primitive.addVertices(5, 8, 7)

        # left bottom
        primitive.addVertices(0, 12, 23)
        primitive.addVertices(23, 12, 0)
        primitive.addVertices(0, 23, 11)
        primitive.addVertices(11, 23, 0)

        # bottom
        primitive.addVertices(0, 1, 13)
        primitive.addVertices(13, 1, 0)
        primitive.addVertices(0, 13, 12)
        primitive.addVertices(12, 13, 0)

        # bottom top left
        primitive.addVertices(11, 23, 22)
        primitive.addVertices(22, 23, 11)
        primitive.addVertices(11, 22, 10)
        primitive.addVertices(10, 22, 11)

        # bottom top right
        primitive.addVertices(3, 15, 14)
        primitive.addVertices(14, 3, 15)
        primitive.addVertices(3, 14, 2)
        primitive.addVertices(2, 14, 3)

        # web left
        primitive.addVertices(10, 22, 21)
        primitive.addVertices(21, 22, 10)
        primitive.addVertices(10, 21, 9)
        primitive.addVertices(9, 21, 10)

        # web right
        primitive.addVertices(3, 16, 15)
        primitive.addVertices(15, 16, 3)
        primitive.addVertices(3, 4, 16)
        primitive.addVertices(16, 4, 3)

        # right bottom
        primitive.addVertices(1, 14, 13)
        primitive.addVertices(13, 14, 1)
        primitive.addVertices(1, 2, 14)
        primitive.addVertices(14, 2, 1)

        # left top
        primitive.addVertices(8, 20, 19)
        primitive.addVertices(19, 20, 8)
        primitive.addVertices(8, 19, 7)
        primitive.addVertices(7, 19, 8)

        # top
        primitive.addVertices(7, 18, 19)
        primitive.addVertices(19, 18, 7)
        primitive.addVertices(7, 18, 6)
        primitive.addVertices(6, 18, 7)

        # top bottom left
        primitive.addVertices(8, 20, 9)
        primitive.addVertices(9, 20, 8)
        primitive.addVertices(20, 9, 21)
        primitive.addVertices(21, 9, 20)

        # top bottom right
        primitive.addVertices(4, 16, 17)
        primitive.addVertices(17, 16, 4)
        primitive.addVertices(5, 16, 17)
        primitive.addVertices(17, 16, 5)
        primitive.addVertices(4, 17, 5)
        primitive.addVertices(5, 17, 4)
        primitive.addVertices(5, 4, 16)
        primitive.addVertices(16, 4, 5)

        # right top
        primitive.addVertices(6, 18, 17)
        primitive.addVertices(17, 18, 6)
        primitive.addVertices(6, 17, 5)
        primitive.addVertices(5, 17, 6)

        # back
        primitive.addVertices(12, 13, 14)
        primitive.addVertices(2+12, 1+12, 0+12)
        primitive.addVertices(0+12, 2+12, 11+12)
        primitive.addVertices(11+12, 2+12, 0+12)
        primitive.addVertices(10+12, 3+12, 4+12)
        primitive.addVertices(4+12, 3+12, 10+12)
        primitive.addVertices(10+12, 4+12, 9+12)
        primitive.addVertices(9+12, 4+12, 10+12)
        primitive.addVertices(6+12, 7+12, 5+12)
        primitive.addVertices(5+12, 7+12, 6+12)
        primitive.addVertices(7+12, 8+12, 5+12)
        primitive.addVertices(5+12, 8+12, 7+12)

        geom = Geom(vertexData)
        geom.addPrimitive(primitive)

        node = GeomNode('I_frame')
        node.addGeom(geom)

        return node
