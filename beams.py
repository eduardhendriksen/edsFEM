import numpy as np


class Beam():
    """
    The Beam class is the representation of physical beams in the structural
    model. It contains all properties and underlying objects needed to model
    a beam.
    """

    def __init__(self, structurenode_1, structurenode_2, system, profile=None,
                 E=None, A=None, I_y=None, I_z=None, G=None, I_t=None,
                 Wy_el=None, Wz_el=None, global_cs=None, a_l=None,
                 sys_eq=False, Wy_pl=None, Wz_pl=None, release=None,
                 segments=12):
        """
        Initiates a Beam object and creates all underlying objects needed to
        construct a Beam object.

        structurenode_1: The starting StructureNode of the Beam
        structurenode_2: The ending StructureNode of the Beam
        system: The System object this Beam is a part of
        E: Young's modulus (N/mm²)
        A: Area (mm²)
        I_y: Moment of Inertia around y-axis (mm^4)
        I_z: Moment of Inertia around z-axis (mm^4)
        G: Shear modulus (N/mm²)
        I_t: torsion property of the section (mm^4)
        Wy_el: Elastic section modulus about the y-axis (mm^3)
        Wz_el: Elastic section modulus about the z-axis (mm^3)
        coord_sys: is the global coordinate system in the form;
                          (np.array([origin]),
                           np.array([[x,x,x],
                                     [y,y,y],
                                     [z,z,z]])
                          in which the vectors x, y and z represent
                          the axes of the sysem
        a_l: the beam's rotation around the longitudinal axis
        Wy_pl: Plastic section modulus about the y-axis
        Wz_pl: Plastic section modulus about the z-axis

        release:  end releases for the element in the form of an numpy array;
                  np.array([x_1, y_1, z_1, r_x_1, r_y_1, r_z_1,
                            x_2, y_2, z_2, r_x_2, r_y_2, r_z_2])
                  by setting the value for the corresponding dof to one, it is
                  released
        segments: The amount of individual elements the beam is seperated in
        """

        self.sys = system

        self.sn_1 = structurenode_1
        self.sn_2 = structurenode_2
        self.point_1 = structurenode_1.node.point
        self.point_2 = structurenode_2.node.point

        self.L = np.sqrt((self.point_2.x - self.point_1.x)**2 +
                         (self.point_2.y - self.point_1.y)**2 +
                         (self.point_2.z - self.point_1.z)**2)

        self.release = release

        self.nodes = []
        self.elements = []

        self.x = None
        self.y = None
        self.z = None

        self.dx = None
        self.dy = None
        self.dz = None

        self.N = None

        self.Mx = None
        self.My = None
        self.Mz = None

        self.Vx = None
        self.Vy = None
        self.Vz = None

        if profile is None:
            self.profile = None
            self.g = 0
            self.designation = 'self assigned'
            self.create_beam(E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                             global_cs=global_cs, a_l=a_l, sys_eq=sys_eq,
                             Wy_pl=Wy_pl, Wz_pl=Wz_pl, release=release,
                             segments=segments)
        else:
            self.profile = profile
            self.g = self.profile.g
            self.designation = self.profile.designation
            if E is None:
                E = profile.E
            if A is None:
                A = profile.A
            if I_y is None:
                I_y = profile.I_y
            if I_z is None:
                I_z = profile.I_z
            if G is None:
                G = profile.G
            if I_t is None:
                I_t = profile.I_t
            if Wy_el is None:
                Wy_el = profile.Wy_el
            if Wz_el is None:
                Wz_el = profile.Wz_el

            self.create_beam(E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                             global_cs=global_cs, a_l=a_l, sys_eq=sys_eq,
                             Wy_pl=Wy_pl, Wz_pl=Wz_pl, release=release,
                             segments=segments)

    def __repr__(self):
        return("\nBeam\n"
               "Point 1: %s\n"
               "Point 2: %s\n"
               "Profile: %s\n" % (self.point_1, self.point_2, self.profile))

    def __str__(self):
        return("\nBeam\n"
               "Point 1: %s\n"
               "Point 2: %s\n"
               "Profile: %s\n" % (self.point_1, self.point_2, self.profile))

    def set_profile(self, profile):
        """
        Changes the beam's profile to the one supplied by changing the element
        properties of the elements in the beam.
        """
        p = profile
        self.profile = p
        self.g = p.g
        self.designation = p.designation
        for el in self.elements:
            el.modify(p.E, p.A, p.I_y, p.I_z, p.G, p.I_t, p.Wy_el, p.Wz_el,
                      Wy_pl=p.Wy_pl, Wz_pl=p.Wz_pl)

    def change_property(self, E=None, A=None, I_y=None, I_z=None, G=None,
                        I_t=None, Wy_el=None, Wz_el=None, Wy_pl=None,
                        Wz_pl=None):
        """
        Change the given properties of the beam by changing the element
        properties of the elements in the beam.
        """
        if E is None:
            if self.profile is not None:
                E = self.profile.E
            else:
                E = self.elements[0].E
        if A is None:
            if self.profile is not None:
                A = self.profile.A
            else:
                A = self.elements[0].A
        if I_y is None:
            if self.profile is not None:
                I_y = self.profile.I_y
            else:
                I_y = self.elements[0].I_y
        if I_z is None:
            if self.profile is not None:
                I_z = self.profile.I_z
            else:
                I_z = self.elements[0].I_z
        if G is None:
            if self.profile is not None:
                G = self.profile.G
            else:
                G = self.elements[0].G
        if I_t is None:
            if self.profile is not None:
                I_t = self.profile.I_t
            else:
                I_t = self.elements[0].I_t
        if Wy_el is None:
            if self.profile is not None:
                Wy_el = self.profile.Wy_el
            else:
                Wy_el = self.elements[0].Wy_el
        if Wz_el is None:
            if self.profile is not None:
                Wz_el = self.profile.Wz_el
            else:
                Wz_el = self.elements[0].Wz_el
        if Wy_pl is None:
            if self.profile is not None:
                try:
                    Wy_pl = self.profile.Wy_pl
                except AttributeError:
                    Wy_pl = None
            else:
                try:
                    Wy_pl = self.elements[0].Wy_pl
                except AttributeError:
                    Wy_pl = None
        if Wz_pl is None:
            if self.profile is not None:
                try:
                    Wz_pl = self.profile.Wz_pl
                except AttributeError:
                    Wz_pl = None
            else:
                try:
                    Wz_pl = self.elements[0].Wz_pl
                except AttributeError:
                    Wz_pl = None

        for el in self.elements:
            el.modify(E, A, I_y, I_z, G, I_t, Wy_el, Wz_el, Wy_pl=Wy_pl,
                      Wz_pl=Wz_pl)

    def create_beam(self, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                    global_cs=None, a_l=None, sys_eq=False, Wy_pl=None,
                    Wz_pl=None, release=None, segments=12):
        """
        Creates a beam with the properties supplied.
        """

        n_1 = self.point_1
        n_2 = self.point_2

        x = n_1.x
        y = n_1.y
        z = n_1.z

        self.x = np.array([x])
        self.y = np.array([y])
        self.z = np.array([z])

        dx = (n_2-n_1).x / segments
        dy = (n_2-n_1).y / segments
        dz = (n_2-n_1).z / segments

        if release is not None:
            self.release = release
            r_1 = np.append(self.release[:6], np.zeros(6))
            r_2 = np.append(np.zeros(6), self.release[6:])
        else:
            self.release = None
            r_1 = None
            r_2 = None

        self.nodes.append(self.sn_1.node)

        c = len(self.sys.nodes) - 1

        for i in range(1, segments):
            if i == 1:

                n = self.sys.add_node((x+i*dx, y+i*dy, z+i*dz))
                c += 1
                el = self.sys.add_element(E, A, I_y, I_z, G, I_t, Wy_el,
                                          Wz_el, self.sn_1.node,
                                          self.sys.nodes[c], sys_eq=sys_eq,
                                          global_cs=global_cs, a_l=a_l,
                                          Wy_pl=Wy_pl, Wz_pl=Wz_pl,
                                          release=r_1, beam=self)

                self.nodes.append(n)
                self.elements.append(el)

                self.x = np.append(self.x, np.array([x+i*dx]))
                self.y = np.append(self.y, np.array([y+i*dy]))
                self.z = np.append(self.z, np.array([z+i*dz]))

            else:

                n = self.sys.add_node((x+i*dx, y+i*dy, z+i*dz))
                c += 1
                el = self.sys.add_element(E, A, I_y, I_z, G, I_t, Wy_el,
                                          Wz_el, self.sys.nodes[c-1],
                                          self.sys.nodes[c], sys_eq=sys_eq,
                                          global_cs=global_cs, a_l=a_l,
                                          Wy_pl=Wy_pl, Wz_pl=Wz_pl,
                                          release=None, beam=self)

                self.nodes.append(n)
                self.elements.append(el)

                self.x = np.append(self.x, np.array([x+i*dx]))
                self.y = np.append(self.y, np.array([y+i*dy]))
                self.z = np.append(self.z, np.array([z+i*dz]))

        el = self.sys.add_element(E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                                  self.sys.nodes[c], self.sn_2.node,
                                  sys_eq=sys_eq, global_cs=global_cs, a_l=a_l,
                                  Wy_pl=Wy_pl, Wz_pl=Wz_pl, release=r_2,
                                  beam=self)

        self.nodes.append(self.sn_2.node)
        self.elements.append(el)

        self.x = np.append(self.x, np.array([self.sn_2.node.point.x]))
        self.y = np.append(self.y, np.array([self.sn_2.node.point.y]))
        self.z = np.append(self.z, np.array([self.sn_2.node.point.z]))

        if self.sys.beams is not None:
            self.sys.beams.append(self)
        else:
            self.sys.beams = [self]

        return self

    def add_node(self, point):
        """
        Adds a node to the beam.
        """

        if self.check_point(point):
            for element in self.elements:
                if point.on_line(element.node_1.point, element.node_2.point):
                    new_el = element.split_by_point(point, self.sys)
                    self.elements.append(new_el)
                    self.sort_elements()
                    return new_el.node_1

    def check_point(self, point):
        """
        Checks if the supplied point is on the beam axis.
        """
        return point.on_line(self.point_1, self.point_2)

    def displacements(self):
        """
        Reads out the displacements in the beams nodes.
        """

        self.dx = None
        self.dy = None
        self.dz = None

        for node in self.nodes:

            u = node.displacement_results

            if type(self.dx) is np.ndarray:
                self.dx = np.append(self.dx, np.array([u[0]]))
            else:
                self.dx = np.array([u[0]])

            if type(self.dy) is np.ndarray:
                self.dy = np.append(self.dy, np.array([u[1]]))
            else:
                self.dy = np.array([u[1]])

            if type(self.dz) is np.ndarray:
                self.dz = np.append(self.dz, np.array([u[2]]))
            else:
                self.dz = np.array([u[2]])

    def moments(self, s=6):
        """
        Reads out the moments in the beams nodes.
        """

        self.Mx = None
        self.My = None
        self.Mz = None

        n = -1

        for element in self.elements:

            n += 1

            x_0 = self.x[n]
            dx = ((self.x[n+1]-self.x[n])/s)
            y_0 = self.y[n]
            dy = ((self.y[n+1]-self.y[n])/s)
            z_0 = self.z[n]
            dz = ((self.z[n+1]-self.z[n])/s)

            d = - 1
            for Mt in element.Mt:
                d += 1
                if type(self.Mx) is np.ndarray:
                    arr = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, Mt])
                    if self.Mx.shape == (4,):
                        self.Mx = np.append(self.Mx.reshape((1, 4)),
                                            arr.reshape((1, 4)), axis=0)
                    else:
                        self.Mx = np.append(self.Mx, arr.reshape((1, 4)),
                                            axis=0)
                else:
                    self.Mx = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, Mt])

            d = - 1
            for My in element.My:
                d += 1
                if type(self.My) is np.ndarray:
                    arr = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, My])
                    if self.My.shape == (4,):
                        self.My = np.append(self.My.reshape((1, 4)),
                                            arr.reshape((1, 4)), axis=0)
                    else:
                        self.My = np.append(self.My, arr.reshape((1, 4)),
                                            axis=0)
                else:
                    self.My = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, My])

            d = - 1
            for Mz in element.Mz:
                d += 1
                if type(self.Mz) is np.ndarray:
                    arr = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, Mz])
                    if self.Mz.shape == (4,):
                        self.Mz = np.append(self.Mz.reshape((1, 4)),
                                            arr.reshape((1, 4)), axis=0)
                    else:
                        self.Mz = np.append(self.Mz, arr.reshape((1, 4)),
                                            axis=0)
                else:
                    self.Mz = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, Mz])

    def normal_forces(self, s=6):
        """
        Reads out the normal forces in the beams nodes.
        """
        self.N = None

        n = -1

        for element in self.elements:

            n += 1

            x_0 = self.x[n]
            dx = ((self.x[n+1]-self.x[n])/s)
            y_0 = self.y[n]
            dy = ((self.y[n+1]-self.y[n])/s)
            z_0 = self.z[n]
            dz = ((self.z[n+1]-self.z[n])/s)

            d = - 1
            for N in element.N:
                d += 1
                if type(self.N) is np.ndarray:
                    arr = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, N])
                    if self.N.shape == (4,):
                        self.N = np.append(self.N.reshape((1, 4)),
                                           arr.reshape((1, 4)), axis=0)
                    else:
                        self.N = np.append(self.N, arr.reshape((1, 4)),
                                           axis=0)
                else:
                    self.N = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, N])

    def shear_forces(self, s=6):
        """
        Reads out the shear forces in the beams nodes.
        """

        self.Vy = None
        self.Vz = None

        n = -1

        for element in self.elements:

            n += 1

            x_0 = self.x[n]
            dx = ((self.x[n+1]-self.x[n])/s)
            y_0 = self.y[n]
            dy = ((self.y[n+1]-self.y[n])/s)
            z_0 = self.z[n]
            dz = ((self.z[n+1]-self.z[n])/s)

            d = - 1
            for Vy in element.Vy:
                d += 1
                if type(self.Vy) is np.ndarray:
                    arr = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, Vy])
                    if self.Vy.shape == (4,):
                        self.Vy = np.append(self.Vy.reshape((1, 4)),
                                            arr.reshape((1, 4)), axis=0)
                    else:
                        self.Vy = np.append(self.Vy, arr.reshape((1, 4)),
                                            axis=0)
                else:
                    self.Vy = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, Vy])

            d = - 1
            for Vz in element.Vz:
                d += 1
                if type(self.Vz) is np.ndarray:
                    arr = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, Vz])
                    if self.Vz.shape == (4,):
                        self.Vz = np.append(self.Vz.reshape((1, 4)),
                                            arr.reshape((1, 4)), axis=0)
                    else:
                        self.Vz = np.append(self.Vz, arr.reshape((1, 4)),
                                            axis=0)
                else:
                    self.Vz = np.array([x_0+d*dx, y_0+d*dy, z_0+d*dz, Vz])

    def sort_elements(self):
        """
        Sorts the elements in the beam from start to end. Also renews the beams
        internal node list.
        """

        if self.elements:
            self.elements = list(dict.fromkeys(self.elements))
            self.elements.sort(key=lambda item:
                               item.node_1.point.on_line(self.point_1,
                                                         self.point_2,
                                                         length=True)[1])

        self.x = np.array([self.point_1.x])
        self.y = np.array([self.point_1.y])
        self.z = np.array([self.point_1.z])
        self.nodes = []

        for el in self.elements:

            if el.node_1 not in self.nodes:
                self.nodes.append(el.node_1)
            if el.node_2 not in self.nodes:
                self.nodes.append(el.node_2)

            self.x = np.append(self.x, np.array(el.node_2.point.x))
            self.y = np.append(self.y, np.array(el.node_2.point.y))
            self.z = np.append(self.z, np.array(el.node_2.point.z))
