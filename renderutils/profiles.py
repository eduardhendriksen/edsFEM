import numpy as np
from panda3d.core import GeomVertexFormat
from panda3d.core import GeomVertexData, GeomVertexWriter, Geom
from panda3d.core import GeomNode, GeomTriangles


class Profile_Renderer():
    """
    Parent class to all profile renderers.
    """

    def __init__(self, base):
        """
        Initiates the Profile_Renderer class
        """

        self.base = base

    def mk_surface(self, p_1, p_2, p_3):
        """
        Method that calculates the surface vertices and normal from the three
        given points.
        """
        u = p_1 - p_2
        v = p_3 - p_2
        n = np.cross(u, v) / np.linalg.norm(np.cross(u, v))

        p = np.append(np.append(np.reshape(p_1, (1, 3)),
                                np.reshape(p_2, (1, 3)), axis=0),
                      np.reshape(p_3, (1, 3)), axis=0)
        return p, n


class I_Profile_Renderer(Profile_Renderer):
    """
    This class is used to render the profile geometry of an I Beam profile.
    """

    def __init__(self, beam, s, base):
        """
        Initiates the I_Profile_Renderer class.
        """
        Profile_Renderer.__init__(self, base)
        self.beam = beam
        self.s = s

    def create_profile(self):
        """
        Creates the necessary rendering geometry to render the 3d rendition of
        an IBeam profile.
        """
        beam = self.beam
        s = self.s
        prof = beam.profile

        c_1 = (0.25882353, 0.80784314, 0.95686275, 1)
        c_2 = (0.19607843, 0.63529412, 0.75686275, 1)

        fmt = GeomVertexFormat.getV3n3c4()
        vertexData = GeomVertexData('I prof', fmt, Geom.UHStatic)

        vertexData.setNumRows(108)

        vertices = GeomVertexWriter(vertexData, 'vertex')
        normals = GeomVertexWriter(vertexData, 'normal')
        colors = GeomVertexWriter(vertexData, 'color')

        T = beam.elements[0].trans_matrix[:3, :3]

        b_0 = np.array([beam.x[0], beam.y[0], beam.z[0]]) / s
        b_l = np.array([beam.x[-1], beam.y[-1], beam.z[-1]]) / s

        points = np.array([[0, prof.b / 2, -prof.h / 2],
                           [0, -prof.b / 2, -prof.h / 2],
                           [0, -prof.b / 2, -prof.h / 2 + prof.t_f],
                           [0, -prof.t_w / 2, -prof.h / 2 + prof.t_f],
                           [0, prof.t_w / 2, -prof.h / 2 + prof.t_f],
                           [0, prof.b / 2, -prof.h / 2 + prof.t_f]
                           ]) / s

        R = np.array([[1, 0, 0],
                      [0, -1, 0],
                      [0, 0, -1]])

        ps = None
        for p in points:
            if ps is None:
                ps = np.reshape(np.dot(R, p), (1, 3))
            else:
                ps = np.append(ps, np.reshape(np.dot(R, p), (1, 3)), axis=0)

        points = np.append(points, ps, axis=0)

        ps = None
        for p in points:
            if ps is None:
                ps = np.reshape(b_0 + np.dot(T.T, p), (1, 3))
            else:
                ps = np.append(ps, np.reshape(b_0 + np.dot(T.T, p), (1, 3)),
                               axis=0)

        pe = None
        for p in points:
            if pe is None:
                pe = np.reshape(b_l + np.dot(T.T, p), (1, 3))
            else:
                pe = np.append(pe, np.reshape(b_l + np.dot(T.T, p), (1, 3)),
                               axis=0)

        points = np.append(ps, pe, axis=0)

        surfs = np.array([[0, 1, 2, c_1], [0, 2, 5, c_1],  # front
                          [4, 3, 10, c_1], [4, 10, 9, c_1],
                          [8, 11, 6, c_1], [8, 6, 7, c_1],
                          [12, 0, 5, c_1], [12, 5, 17, c_1],  # t./b. sides
                          [1, 13, 14, c_1], [1, 14, 2, c_1],
                          [11, 23, 18, c_1], [11, 18, 6, c_1],
                          [20, 8, 7, c_1], [20, 7, 19, c_1],
                          [16, 4, 9, c_2], [16, 9, 21, c_2],  # Mid sides
                          [3, 15, 22, c_2], [3, 22, 10, c_2],
                          [13, 12, 17, c_1], [13, 17, 14, c_1],  # Back
                          [15, 16, 21, c_1], [15, 21, 22, c_1],
                          [23, 20, 19, c_1], [23, 19, 18, c_1],
                          [7, 6, 18, c_1], [7, 18, 19, c_1],  # T. t. fl.
                          [9, 8, 20, c_1], [9, 20, 21, c_1],  # In. t. fl. l
                          [11, 10, 22, c_1], [11, 22, 23, c_1],  # In. t. fl. r
                          [5, 4, 16, c_1], [5, 16, 17, c_1],  # In. b. fl. l.
                          [3, 2, 14, c_1], [3, 14, 15, c_1],  # In. b. fl. l.
                          [1, 0, 12, c_1], [1, 12, 13, c_1],  # Bot. bot. fl.
                          ])

        prim = GeomTriangles(Geom.UHStatic)

        vtx_count = -1

        for surf in surfs:
            p, n = self.mk_surface(points[surf[0]],
                                   points[surf[1]],
                                   points[surf[2]])

            for i in range(p.shape[0]):
                vertices.addData3f(*p[i, :])
                normals.addData3f(*(-n))
                colors.addData4f(*surf[3])
                vtx_count += 1

            prim.add_vertices(vtx_count - 2,
                              vtx_count - 1,
                              vtx_count)

        geom = Geom(vertexData)
        geom.addPrimitive(prim)

        node = GeomNode('I_frame')
        node.addGeom(geom)
        return node


class CHS_Renderer(Profile_Renderer):
    """
    This class is used to render the profile geometry of a Circular Hollow
    Section or pipe.
    """

    def __init__(self, beam, s, base):
        """
        Initiates the CHS_Renderer class.
        """
        Profile_Renderer.__init__(self, base)
        self.beam = beam
        self.s = s

    def create_profile(self):
        """
        Creates the necessary rendering geometry to render the 3d rendition of
        a CHS.
        """

        beam = self.beam
        s = self.s
        prof = beam.profile

        color = (0.25882353, 0.80784314, 0.95686275, 1)

        fmt = GeomVertexFormat.getV3n3c4()
        vertexData = GeomVertexData('something', fmt, Geom.UHStatic)

        vertexData.setNumRows(2304)

        vertices = GeomVertexWriter(vertexData, 'vertex')
        normals = GeomVertexWriter(vertexData, 'normal')
        colors = GeomVertexWriter(vertexData, 'color')

        T = beam.elements[0].trans_matrix[:3, :3]

        b_0 = np.array([beam.x[0], beam.y[0], beam.z[0]]) / s
        b_l = np.array([beam.x[-1], beam.y[-1], beam.z[-1]]) / s

        s_ang = 49

        wall = []

        for in_out in ((prof.D / 2) / s, (prof.D / 2 - prof.t) / s):
            for b in [b_0, b_l]:
                layer = []
                for j in range(s_ang):
                    ang = j * (2 * np.pi / (s_ang-1))

                    p = b + np.dot(T.T, np.array([0,
                                                  in_out * np.sin(ang),
                                                  in_out * np.cos(ang)]))

                    layer.append(p)
                if layer:
                    wall.append(layer)

        prim = GeomTriangles(Geom.UHStatic)

        vtx_count = -1

        for io in (0, 1):
            for l in range(2):
                for a in range(s_ang-1):
                    if io == 0:
                        p, n = self.mk_surface(wall[l][a],
                                               wall[l+1][a+1],
                                               wall[l][a+1])
                    else:
                        p, n = self.mk_surface(wall[l+io][a],
                                               wall[l+io][a+1],
                                               wall[l+1+io][a+1])

                    for i in range(p.shape[0]):
                        vertices.addData3f(*p[i, :])
                        normals.addData3f(*(-n))
                        colors.addData4f(*color)
                        vtx_count += 1

                    prim.add_vertices(vtx_count - 2,
                                      vtx_count - 1,
                                      vtx_count)

                    if io == 0:
                        p, n = self.mk_surface(wall[l+io][a],
                                               wall[l+1+io][a],
                                               wall[l+1+io][a+1])
                    else:
                        p, n = self.mk_surface(wall[l+io][a],
                                               wall[l+1+io][a+1],
                                               wall[l+1+io][a])

                    for i in range(p.shape[0]):
                        vertices.addData3f(*p[i, :])
                        normals.addData3f(*(-n))
                        colors.addData4f(*color)
                        vtx_count += 1

                    prim.add_vertices(vtx_count - 2,
                                      vtx_count - 1,
                                      vtx_count)

        for io in (0, 2):
            for a in range(s_ang-1):
                if io == 0:
                    p, n = self.mk_surface(wall[0][a],
                                           wall[-2][a+1],
                                           wall[-2][a])
                else:
                    p, n = self.mk_surface(wall[1][a],
                                           wall[-1][a],
                                           wall[-1][a+1])
                for i in range(p.shape[0]):
                    vertices.addData3f(*p[i, :])
                    normals.addData3f(*(-n))
                    colors.addData4f(*color)
                    vtx_count += 1
                prim.add_vertices(vtx_count - 2,
                                  vtx_count - 1,
                                  vtx_count)
                if io == 0:
                    p, n = self.mk_surface(wall[0][a],
                                           wall[0][a+1],
                                           wall[-2][a+1])
                else:
                    p, n = self.mk_surface(wall[1][a],
                                           wall[-1][a+1],
                                           wall[1][a+1])
                for i in range(p.shape[0]):
                    vertices.addData3f(*p[i, :])
                    normals.addData3f(*(-n))
                    colors.addData4f(*color)
                    vtx_count += 1
                prim.add_vertices(vtx_count - 2,
                                  vtx_count - 1,
                                  vtx_count)

        geom = Geom(vertexData)
        geom.addPrimitive(prim)

        node = GeomNode('CHS')
        node.addGeom(geom)
        return node
