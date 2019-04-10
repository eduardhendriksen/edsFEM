import numpy as np
from scipy.linalg import expm


class BeamElement:
    """
    The BeamElement class is the elemental representation of portion of a
    structural beam.
    """

    def __init__(self, ID, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                 node_1, node_2, global_cs=None, a_l=None, sys_eq=False,
                 Wy_pl=None, Wz_pl=None, release=None, beam=None):
        """
        ID: integer representing the elements ID
        E: Young's modulus
        A: Area
        I_y: Moment of Inertia around y-axis
        I_z: Moment of Inertia around z-axis
        G: Shear modulus
        I_t: torsion property of the section
        Wy_el: Elastic section modulus about the y-axis
        Wz_el: Elastic section modulus about the z-axis
        node_1: node object
        node_2: node object
        coord_sys: is the global coordinate system in the form;
                          (np.array([origin]),
                           np.array([[x,x,x],
                                     [y,y,y],
                                     [z,z,z]])
                          in which the vectors x, y and z represent
                          the axes of the sysem
        Wy_pl: Plastic section modulus about the y-axis
        Wz_pl: Plastic section modulus about the z-axis
        a_l: the beam's rotation around the longitudinal axis in degrees
        release:  end releases for the element in the form of an numpy array;
                  np.array([x_1, y_1, z_1, r_x_1, r_y_1, r_z_1,
                            x_2, y_2, z_2, r_x_2, r_y_2, r_z_2])
                  by setting the value for the corresponding dof to one, it is
                  released
        """
        self.ID = ID
        self.E = E
        self.A = A
        self.I_z = I_z
        self.I_y = I_y
        self.EA = self.E*self.A
        self.EI_y = self.E*self.I_y
        self.EI_z = self.E*self.I_z
        self.G = G
        self.I_t = I_t
        self.GI_t = self.G*self.I_t
        self.Wy_el = Wy_el
        self.Wz_el = Wz_el
        self.Wy_pl = Wy_pl
        self.Wz_pl = Wz_pl

        self.node_1 = node_1

        if self.node_1.elements is not None:
            if self not in self.node_1.elements:
                self.node_1.elements.append(self)
        else:
            self.node_1.elements = [self]

        self.node_2 = node_2

        if self.node_2.elements is not None:
            if self not in self.node_2.elements:
                self.node_2.elements.append(self)
        else:
            self.node_2.elements = [self]

        self.point_1 = None
        self.point_2 = None
        self.L = None

        self.length()

        self.global_cs = global_cs
        self.local_cs = None
        self.sys_eq = sys_eq

        if a_l is not None:
            self.a_l = (a_l/180)*np.pi
        else:
            self.a_l = a_l

        self.release = release

        self.beam = beam

        self.coord_sys()
        self.trans_matrix = None
        self.trans_m()
        self.loc_stiffness_matrix = None
        self.glob_stiffness_matrix = None

        self.displacements = None
        self.forces = None

    def length(self):
        """
        Calculates the element length based on its start and end point
        """

        self.point_1 = self.node_1.point
        self.point_2 = self.node_2.point

        self.L = np.sqrt((self.point_2.x - self.point_1.x)**2 +
                         (self.point_2.y - self.point_1.y)**2 +
                         (self.point_2.z - self.point_1.z)**2)

    def coord_sys(self):
        """
        The local coordinate system always has its x-axis along the
        longitudinal axis of the beam element.
        If a_l is not specified, the z-axis is upward.

               z ↺ ↖
           ↙ |__ y ↺
             /
            x  ↗
           ↺

        Method to transform coordinate systems is taken from;
        'https://math.stackexchange.com/questions/180418/calculate-rotation-
         matrix-to-align-vector-a-to-vector-b-in-3d'
        """
        if self.global_cs is None:

            self.global_cs = (np.array([0, 0, 0]),
                              np.array([[1, 0, 0],
                                        [0, 1, 0],
                                        [0, 0, 1]]))

        else:

            self.global_cs = (self.global_cs[0],
                              np.array([(self.global_cs[1][0] /
                                         np.linalg.norm(self.global_cs[1][0])),
                                        (self.global_cs[1][1] /
                                         np.linalg.norm(self.global_cs[1][1])),
                                        (self.global_cs[1][2] /
                                         np.linalg.norm(self.global_cs[1][2]))
                                        ]))

        if self.sys_eq is True:

            self.local_cs = self.global_cs

        else:

            origin = np.array([self.point_1.x, self.point_1.y, self.point_1.z])

            x_x = self.point_2.x - self.point_1.x
            x_y = self.point_2.y - self.point_1.y
            x_z = self.point_2.z - self.point_1.z

            x_axis = np.array([x_x, x_y, x_z])
            x_axis = x_axis / np.linalg.norm(x_axis)

            if x_axis[2] == 0:

                z_axis = np.array([0, 0, 1])

            elif x_axis[0] == x_axis[1] == 0:

                z_axis = np.array([-1, 0, 0])

            else:

                z_dir = (x_axis[0]**2 + x_axis[1]**2) / x_axis[2]

                z_axis = (np.array([-x_axis[0], -x_axis[1], z_dir]) /
                          np.linalg.norm(np.array([-x_axis[0], -x_axis[1],
                                                   z_dir])))

            if np.isclose(np.dot(x_axis, z_axis), 0):

                y_axis = - np.cross(x_axis, z_axis)

            else:

                raise ValueError("Can't create an orthogonal coordinate system"
                                 + " at element " + str(self.ID))

            if self.a_l is not None:

                y_axis = np.dot(self.rot_matrix(self.a_l, x_axis), y_axis)
                z_axis = np.dot(self.rot_matrix(self.a_l, x_axis), z_axis)

            if (np.isclose(np.dot(x_axis, z_axis), 0) and
                np.isclose(np.dot(x_axis, y_axis), 0) and
                    np.isclose(np.dot(y_axis, z_axis), 0)):

                self.local_cs = (origin,
                                 np.append(np.append(x_axis.reshape((1, 3)),
                                                     y_axis.reshape((1, 3)),
                                                     axis=0),
                                           z_axis.reshape((1, 3)), axis=0))
            else:

                raise ValueError('Local coordinate system not orthognal for '
                                 'element ' + str(self.ID))

    def trans_m(self):
        """
        Create the transformation matrix between the local and global
        coordinate system using direction cosines (which are the normalized
        local axis vectors calculated in the coord_sys() function)

        T.T k_loc T = k_glob
        T   u_glob  = u_loc
        T.T F_loc   = F_glob
        T.T loc_axis = glob_axis

        """

        if self.local_cs and self.global_cs:

            dcs = self.local_cs[1]

            row_1 = np.append(dcs, np.zeros((3, 9)), axis=1)
            row_2 = np.append(np.append(np.zeros((3, 3)), dcs, axis=1),
                              np.zeros((3, 6)), axis=1)
            row_3 = np.append(np.append(np.zeros((3, 6)), dcs, axis=1),
                              np.zeros((3, 3)), axis=1)
            row_4 = np.append(np.zeros((3, 9)), dcs, axis=1)

            self.trans_matrix = np.append(np.append(np.append(row_1, row_2,
                                                              axis=0),
                                                    row_3, axis=0),
                                          row_4, axis=0)

        else:

            raise ValueError('No local or global cs specified at element ' +
                             str(self.ID))

    def rot_matrix(self, angle, axis):
        """
        Returns a 3d rotation matrix based on an axis np.array([x,y,z])
        and the angle (in radians) of the counterclockwise rotation around
        said axis
        """
        return expm(np.cross(np.eye(3), np.asarray(axis) /
                             np.linalg.norm(np.asarray(axis))*angle))

    def split_by_point(self, point, sys):
        """
        Splits the element into two elements by a given point.
        """

        if point.on_line(self.node_1.point, self.node_2.point):

            if type(self) is EBElement:

                nn = sys.add_node(point)

                if self.release is not None:

                    r_1 = np.append(self.release[:6], np.zeros(6))
                    r_2 = np.append(np.zeros(6), self.release[6:])

                else:

                    r_1 = None
                    r_2 = None

                ne = sys.add_element(self.E, self.A, self.I_y, self.I_z,
                                     self.G, self.I_t, self.Wy_el,
                                     self.Wz_el, nn, self.node_2,
                                     global_cs=self.global_cs, a_l=self.a_l,
                                     sys_eq=self.sys_eq, Wy_pl=self.Wy_pl,
                                     Wz_pl=self.Wz_pl, release=r_2,
                                     beam=self.beam)

                self.node_2.elements.remove(self)

                self.node_2 = nn

                if self.node_2.elements is not None:

                    if self not in self.node_2.elements:

                        self.node_2.elements.append(self)

                else:

                    self.node_2.elements = [self]

                self.length()
                self.coord_sys()
                self.trans_m()
                self.local_stiffness_matrix()
                self.global_stiffness_matrix()
                self.releases(release=r_1)

                self.beam.elements.append(ne)
                self.beam.sort_elements()

                return ne

            else:

                raise NotImplementedError('Only possible with EB Elements')
        else:

            print()
            print('point: ' + str(point.v))
            print()
            print('el_p1: ' + str(self.node_1.point.v))
            print('el_p2: ' + str(self.node_2.point.v))
            print()
            raise ValueError("Point not in element")

    def internal_forces(self, s=6):
        """
        Uses the element local stiffness matrix, element transformation matrix
        and element global displacements to determine the internal element
        forces in the local coordinate system.

        Also assigns global force data to support nodes.
        """

        u = self.node_1.displacement_results
        u = np.append(u, self.node_2.displacement_results)
        self.displacements = u

        if any(k is not None and k != 0 for n in [self.node_1, self.node_2]
               for k in [n.support, n.Fx, n.Fy, n.Fz, n.Tx, n.Ty, n.Tz]):

            f_glob = np.dot(self.glob_stiffness_matrix, self.displacements)

            if self.node_1.support is not None:

                self.node_1.force_results = f_glob[:6]

                if [k for k in [self.node_1.Fx,
                                self.node_1.Fy,
                                self.node_1.Fz,
                                self.node_1.Tx,
                                self.node_1.Ty,
                                self.node_1.Tz] if k != 0]:

                    self.node_1.force_results += np.array([self.node_1.Fx,
                                                           self.node_1.Fy,
                                                           self.node_1.Fz,
                                                           self.node_1.Tx,
                                                           self.node_1.Ty,
                                                           self.node_1.Tz])
            if self.node_2.support is not None:

                self.node_2.force_results = f_glob[6:]

            elif [k for k in [self.node_2.Fx,
                              self.node_2.Fy,
                              self.node_2.Fz,
                              self.node_2.Tx,
                              self.node_2.Ty,
                              self.node_2.Tz] if k != 0]:

                self.node_2.force_results = np.array([self.node_2.Fx,
                                                      self.node_2.Fy,
                                                      self.node_2.Fz,
                                                      self.node_2.Tx,
                                                      self.node_2.Ty,
                                                      self.node_2.Tz])

        u = np.dot(self.trans_matrix, u)
        K = self.loc_stiffness_matrix
        f = np.dot(K, u)
        self.forces = f

        self.normal_force(s=s)
        self.shear_force(s=s)
        self.moment(s=s)
        self.torsion(s=s)

        self.S = []

        T = self.trans_matrix[:3, :3]

        for i in range(s):

            if type(self.S) is np.ndarray:

                self.S = np.append(self.S,
                                   (self.S[-1, :] +
                                    np.dot(T.T, np.array([self.L/s, 0, 0]))
                                    ).reshape((1, 3)), axis=0)

            else:

                self.S = (np.array([0, 0, 0]) + self.local_cs[0]
                          ).reshape((1, 3))

        self.max_stresses()

    def normal_force(self, s=6):
        """
        Determines the axial normal force in the element.

        N positive  ← ███ →

        """
        n_1 = self.forces[0]
        n_2 = self.forces[6]

        self.N = []

        if np.isclose(np.abs(n_1), np.abs(n_2)):

            self.N = np.full((s+1,), (n_2-n_1) / 2)

        else:

            for i in range(s+1):

                if type(self.N) is np.ndarray:

                    self.N = np.append(self.N, [self.N[-1] + (n_2 + n_1) / s])

                else:

                    self.N = np.array([-n_1])

    def shear_force(self, s=6):
        """
        Determines the shear force in the element in the two shear directions.

        Vy positive  ↑ ███ ↓

        Vz positive  ↓ ███ ↑

        """
        v_y_1 = self.forces[1]
        v_y_2 = self.forces[7]

        self.Vy = []

        if np.isclose(np.abs(v_y_1), np.abs(v_y_2)):

            self.Vy = np.full((s+1,), (v_y_1 - v_y_2) / 2)

        else:

            for i in range(s+1):

                if type(self.Vy) is np.ndarray:

                    self.Vy = np.append(self.Vy, [self.Vy[-1] +
                                                  -(v_y_2 + v_y_1) / s])

                else:

                    self.Vy = np.array([v_y_1])

        v_z_1 = self.forces[2]
        v_z_2 = self.forces[8]

        self.Vz = []

        if np.isclose(np.abs(v_z_1), np.abs(v_z_2)):

            self.Vz = np.full((s+1,), (v_z_2 - v_z_1) / 2)

        else:

            for i in range(s+1):

                if type(self.Vz) is np.ndarray:

                    self.Vz = np.append(self.Vz, [self.Vz[-1] +
                                                  (v_z_2 + v_z_1) / s])

                else:

                    self.Vz = np.array([-v_z_1])

    def moment(self, s=6):
        """
        Determines the bending moment in the element in the two bending
        directions.

        My positive ↳ ███ ↲

        Mz positive ↱ ███ ↰

        """
        m_y_1 = self.forces[4]
        m_y_2 = self.forces[10]

        self.My = []

        if np.isclose(np.abs(m_y_1), np.abs(m_y_2)):

            self.My = np.full((s+1,), (m_y_2 - m_y_1) / 2)

        else:

            for i in range(s+1):

                if type(self.My) is np.ndarray:

                    self.My = np.append(self.My, [self.My[-1] +
                                                  (m_y_2 + m_y_1) / s])

                else:

                    self.My = np.array([-m_y_1])

        m_z_1 = self.forces[5]
        m_z_2 = self.forces[11]

        self.Mz = []

        if np.isclose(np.abs(m_z_1), np.abs(m_z_2)):

            self.Mz = np.full((s+1,), (m_z_1 - m_z_2) / 2)

        else:

            for i in range(s+1):

                if type(self.Mz) is np.ndarray:

                    self.Mz = np.append(self.Mz, [self.Mz[-1] +
                                                  -(m_z_2 + m_z_1) / s])

                else:

                    self.Mz = np.array([m_z_1])

    def torsion(self, s=6):
        """
        Determines the torsion moment in the element in the element
        longitudinal axis.

        positive torsion   --- x ↺

        """

        m_t_1 = self.forces[3]
        m_t_2 = self.forces[9]

        self.Mt = []

        if np.isclose(np.abs(m_t_1), np.abs(m_t_2)):

            self.Mt = np.full((s+1,), (m_t_2 - m_t_1) / 2)

        else:

            for i in range(s+1):

                if type(self.Mt) is np.ndarray:

                    self.Mt = np.append(self.Mt, [self.Mt[-1] +
                                                  (m_t_2 + m_t_1) / s])

                else:

                    self.Mt = np.array([-m_t_1])

    def max_stresses(self):
        """
        Determines the maximal shear and normal stresses due to the
        internal forces in the element, not taking into account torsion.
        """

        self.sigma = []
        self.tau = []
        S_z = 0

        if self.beam:

            if self.beam.profile:

                try:

                    p = self.beam.profile
                    S_z = (0.5 * (p.h-2*p.t_f) * p.t_w) * 0.25 * (p.h-2*p.t_f)
                    S_z += p.b * p.t_f * (0.5 * p.h - 0.5 * p.t_f)
                    S_y = ((p.b-0.5*p.t_w) * p.t_f) * 0.25 * (p.b)
                    S_y += (p.h-2*p.t_f)*0.5*p.t_w*0.25*p.t_w

                except AttributeError:

                    pass

        for i in range(self.S.shape[0]):

            sigma_normal = np.max(np.abs(self.N[i] / self.A))
            m_y = np.max(np.abs(self.My[i])) / self.Wy_el
            m_z = np.max(np.abs(self.Mz[i])) / self.Wz_el
            sigma_bending = np.sum([m_y, m_z])

            if type(self.sigma) is np.ndarray:

                self.sigma = np.append(self.sigma,
                                       np.array([sigma_normal +
                                                 sigma_bending])
                                       )
            else:

                self.sigma = np.array([sigma_normal + sigma_bending])

            if S_z > 0:

                tau_mz = (self.Vz[i] * S_z) / (p.I_y * p.t_w)
                tau_my = (self.Vy[i] * S_y) / (p.I_z * p.t_w)

                tau = np.abs(tau_mz) + np.abs(tau_my)

                if type(self.tau) is np.ndarray:

                    self.tau = np.append(self.tau, np.array(tau))

                else:

                    self.tau = np.array([tau])

    def releases(self, release=None):
        """
        Employs static condensation to implement member end releases.

        release:  end releases for the element in the form of an numpy array;
                  np.array([x_1, y_1, z_1, r_x_1, r_y_1, r_z_1,
                            x_2, y_2, z_2, r_x_2, r_y_2, r_z_2])
                  by setting the value for the corresponding dof to one, it is
                  released
        """

        if release is None:

            if self.release is None:

                return

            elif np.all(np.asarray(self.release) == np.zeros(12)):

                return

            else:

                release = np.array(self.release)

        elif np.all(np.asarray(release) == np.zeros(12)):

            return

        self.local_stiffness_matrix()
        k = self.loc_stiffness_matrix

        n_1 = self.node_1
        n_2 = self.node_2

        f = np.array([n_1.Fx, n_1.Fy, n_1.Fz, n_1.Tx, n_1.Ty, n_1.Tz,
                      n_2.Fx, n_2.Fy, n_2.Fz, n_2.Tx, n_2.Ty, n_2.Tz])

        release = np.asarray(release)
        cd = np.where(release == 1)[0]

        aindx = np.arange(k.shape[0])
        aindx = np.delete(aindx, cd, 0)

        bindx = cd
        Kaa = np.mat(k[np.ix_(aindx, aindx)])
        Kab = np.mat(k[np.ix_(aindx, bindx)])
        Kbb = np.mat(k[np.ix_(bindx, bindx)])

        fa = np.mat(f[aindx])
        fb = np.mat(f[bindx])

        K1 = Kaa-Kab*Kbb.I*Kab.T

        k_ = np.zeros(k.shape)

        for n, i in zip(aindx, np.arange(K1.shape[0])):

            for l, m in zip(aindx, np.arange(K1.shape[0])):

                k_[n, l] = K1[i, m]

        self.loc_stiffness_matrix = k_

        f1 = fa.T-np.dot(Kab, np.dot(Kbb.I, fb.T))
        f_ = np.zeros(k.shape[0])

        for n, i in zip(aindx, np.arange(f1.shape[0])):

            f_[n] = f1[i]

        (n_1.Fx, n_1.Fy, n_1.Fz, n_1.Tx, n_1.Ty, n_1.Tz,
         n_2.Fx, n_2.Fy, n_2.Fz, n_2.Tx, n_2.Ty, n_2.Tz) = f_

        self.loc_stiffness_matrix = k_

        self.global_stiffness_matrix()

    def modify(self, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
               Wy_pl=None, Wz_pl=None):
        """
        Changes the element properties to the ones supplied. Recalculates the
        local and global stiffness matrices and reapplies releases.
        """
        self.E = E
        self.A = A
        self.E = E
        self.A = A
        self.I_z = I_z
        self.I_y = I_y
        self.EA = self.E*self.A
        self.EI_y = self.E*self.I_y
        self.EI_z = self.E*self.I_z
        self.G = G
        self.I_t = I_t
        self.GI_t = self.G*self.I_t
        self.Wy_el = Wy_el
        self.Wz_el = Wz_el
        self.Wy_pl = Wy_pl
        self.Wz_pl = Wz_pl

        self.local_stiffness_matrix()
        self.global_stiffness_matrix()
        self.releases()


class EBElement(BeamElement):
    """
    The EBElement is Euler-Bernoulli representation of the beam element.
    """

    def __init__(self, ID, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                 node_1, node_2, global_cs=None, a_l=None, sys_eq=False,
                 Wy_pl=None, Wz_pl=None, release=None, beam=None):
        """
        ID: integer representing the elements ID
        E: Young's modulus
        A: Area
        I_y: Moment of Inertia around y-axis
        I_z: Moment of Inertia around z-axis
        G: Shear modulus
        I_t: torsion property of the section
        Wy_el: Elastic section modulus about the y-axis
        Wz_el: Elastic section modulus about the z-axis
        node_1: node object
        node_2: node object
        coord_sys: is the global coordinate system in the form;
                          (np.array([origin]),
                           np.array([[x,x,x],
                                     [y,y,y],
                                     [z,z,z]])
                          in which the vectors x, y and z represent
                          the axes of the sysem
        Wy_pl: Plastic section modulus about the y-axis
        Wz_pl: Plastic section modulus about the z-axis
        a_l: the beam's rotation around the longitudinal axis
        release:  end releases for the element in the form of an numpy array;
                  np.array([x_1, y_1, z_1, r_x_1, r_y_1, r_z_1,
                            x_2, y_2, z_2, r_x_2, r_y_2, r_z_2])
                  by setting the value for the corresponding dof to one, it is
                  released
        """
        BeamElement.__init__(self, ID, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                             node_1, node_2, global_cs=global_cs,
                             a_l=a_l, sys_eq=sys_eq, Wy_pl=Wy_pl,
                             Wz_pl=Wz_pl, release=release, beam=beam)

        self.local_stiffness_matrix()
        self.global_stiffness_matrix()

    def local_stiffness_matrix(self):
        """
        This functions assembles the local stiffness matrix according to the
        EB_K_matrix.png file in the package.
        http://www.serendi-cdi.org/serendipedia/index.php?title=Stiffness

        With the change that the torsion terms depend on I_t, the empirical
        torsion constant for the beam, and not J, it's polar moment of inertia.
        This is due to the fact that structural steel beams rarely display
        the torsion capacity associated with the polar moment of inertia.
        The torsion resistance is usually much lower.
        """

        m = np.zeros((12, 12))

        # X_1
        m[0, 0] = self.EA / self.L
        m[0, 6] = -self.EA / self.L
        # Y_1
        m[1, 1] = 12 * self.EI_z / self.L**3
        m[1, 5] = 6 * self.EI_z / self.L**2
        m[1, 7] = -12 * self.EI_z / self.L**3
        m[1, 11] = 6 * self.EI_z / self.L**2
        # Z_1
        m[2, 2] = 12 * self.EI_y / self.L**3
        m[2, 4] = -6 * self.EI_y / self.L**2
        m[2, 8] = -12 * self.EI_y / self.L**3
        m[2, 10] = -6 * self.EI_y / self.L**2
        # t_X_1
        m[3, 3] = self.GI_t / self.L
        m[3, 9] = -self.GI_t / self.L
        # t_Y_1
        m[4, 2] = -6 * self.EI_y / self.L**2
        m[4, 4] = 4 * self.EI_y / self.L
        m[4, 8] = 6 * self.EI_y / self.L**2
        m[4, 10] = 2 * self.EI_y / self.L
        # t_Z_1
        m[5, 1] = 6 * self.EI_z / self.L**2
        m[5, 5] = 4 * self.EI_z / self.L
        m[5, 7] = -6 * self.EI_z / self.L**2
        m[5, 11] = 2 * self.EI_z / self.L
        # X_2
        m[6, 0] = -self.EA / self.L
        m[6, 6] = self.EA / self.L
        # Y_2
        m[7, 1] = -12 * self.EI_z / self.L**3
        m[7, 5] = -6 * self.EI_z / self.L**2
        m[7, 7] = 12 * self.EI_z / self.L**3
        m[7, 11] = -6 * self.EI_z / self.L**2
        # Z_2
        m[8, 2] = -12 * self.EI_y / self.L**3
        m[8, 4] = 6 * self.EI_y / self.L**2
        m[8, 8] = 12 * self.EI_y / self.L**3
        m[8, 10] = 6 * self.EI_y / self.L**2
        # t_X_2
        m[9, 3] = -self.GI_t / self.L
        m[9, 9] = self.GI_t / self.L
        # t_Y_2
        m[10, 2] = -6 * self.EI_y / self.L**2
        m[10, 4] = 2 * self.EI_y / self.L
        m[10, 8] = 6 * self.EI_y / self.L**2
        m[10, 10] = 4 * self.EI_y / self.L
        # t_Z_2
        m[11, 1] = 6 * self.EI_z / self.L**2
        m[11, 5] = 2 * self.EI_z / self.L
        m[11, 7] = -6 * self.EI_z / self.L**2
        m[11, 11] = 4 * self.EI_z / self.L

        if np.allclose(m, m.T):
            self.loc_stiffness_matrix = m
        else:
            raise ValueError('The local stiffness matrix for element ' +
                             str(self.ID) + ' is not symmetric')

    def global_stiffness_matrix(self):
        """
        This method transforms the local stiffness matrix into its global
        counterpart
        """
        if self.loc_stiffness_matrix is not None:

            if self.trans_matrix is not None:

                T = self.trans_matrix
                k = self.loc_stiffness_matrix
                k_ = np.dot(T.T, np.dot(k, T))

                if np.allclose((k_.transpose()), k_):

                    self.glob_stiffness_matrix = k_

                else:

                    print()
                    print()
                    print('Biggest difference between k.T and k')
                    print(np.max(k_.transpose() - k_))
                    print()
                    print()

                    raise ValueError('The global stiffness matrix for element '
                                     + str(self.ID) + ' is not symmetric')
            else:

                raise ValueError('Define the transformation matrix before ' +
                                 'transforming the stiffness matrix of ' +
                                 'element ' + str(self.ID))

        else:

            raise ValueError('The local stiffness matrix for element ' +
                             str(self.ID) + ' is not present')
