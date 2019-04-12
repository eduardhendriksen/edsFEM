import numpy as np
from edsFEM.points import Point
from edsFEM.rectangles import Rectangle


class PointLoad():
    """
    The PointLoad class contains the node the loads are applied to and
    methods to apply and remove the loads.
    """

    def __init__(self, node, name=None, Fx=0, Fy=0, Fz=0, Tx=0,
                 Ty=0, Tz=0):
        """
        Initializes the PointLoad class

        node: Node class the point load is applied to
        name: Name of the load
        Fx: Load in global x-direction (N)
        Fy: Load in global y-direction (N)
        Fz: Load in global z-direction (N)
        Tx: Moment about the global x-axis (Nmm)
        Ty: Moment about the global y-axis (Nmm)
        Tz: Moment about the global z-axis (Nmm)
        """
        self.node = node
        self.name = name
        self.Fx = Fx
        self.Fy = Fy
        self.Fz = Fz
        self.Tx = Tx
        self.Ty = Ty
        self.Tz = Tz

    def __repr__(self):
        return("\nPoint Load\n"
               "Node: %s\n"
               "Name: %s\n"
               "Fx: %s\n"
               "Fy: %s\n"
               "Fz: %s\n"
               "Tx: %s\n"
               "Ty: %s\n"
               "Tz: %s\n" % (self.node,
                             self.name,
                             round(self.Fx, 2),
                             round(self.Fy, 2),
                             round(self.Fz, 2),
                             round(self.Tx, 2),
                             round(self.Ty, 2),
                             round(self.Tz, 2)))

    def __str__(self):
        return("\nPoint Load\n"
               "Node: %s\n"
               "Name: %s\n"
               "Fx: %s\n"
               "Fy: %s\n"
               "Fz: %s\n"
               "Tx: %s\n"
               "Ty: %s\n"
               "Tz: %s\n" % (self.node,
                             self.name,
                             round(self.Fx, 2),
                             round(self.Fy, 2),
                             round(self.Fz, 2),
                             round(self.Tx, 2),
                             round(self.Ty, 2),
                             round(self.Tz, 2)))

    def apply(self):
        """
        This method applies the point load to the node.
        """
        self.node.Fx += self.Fx
        self.node.Fy += self.Fy
        self.node.Fz += self.Fz
        self.node.Tx += self.Tx
        self.node.Ty += self.Ty
        self.node.Tz += self.Tz

    def remove(self):
        """
        This method removes the point load from the node.
        """
        self.node.Fx -= self.Fx
        self.node.Fy -= self.Fy
        self.node.Fz -= self.Fz
        self.node.Tx -= self.Tx
        self.node.Ty -= self.Ty
        self.node.Tz -= self.Tz


class QLoad():
    """
    The QLoad class contains all the properties and methods to apply
    distributed loads to the structure.
    """

    def __init__(self, start_val, start_point, end_point, name=None,
                 end_val=None, direction=[0, 0, -1], beam=None):
        """
        start_val : in kN/m
        start_point : point object
        end_point : point object
        end_val : optional (linearly changing load) (kN / m )
        direction : normalized vector of q_load direction, standard in
                    -z direction (global coordinate system)
        beam: optional beam the load is constrained to
        """
        if name is not None:
            self.name = name
        else:
            self.name = 'q_load'
        self.start_val = start_val
        if end_val is not None:
            self.end_val = end_val
        else:
            self.end_val = start_val

        self.start_point = start_point
        self.end_point = end_point
        self.direction = np.array(direction)
        self.beam = beam
        self.L = np.sqrt(np.square(np.sum([self.start_point.x-self.end_point.x,
                                           self.start_point.y-self.end_point.y,
                                           self.start_point.z-self.end_point.z]
                                          )))
        self.total_load = (self.start_val * self.L +
                           (self.end_val - self.start_val) * self.L / 2)
        self.s_load = self.start_val * self.direction
        self.e_load = self.end_val * self.direction
        self.d_load = (self.e_load - self.s_load) / self.L
        self.element_list = []
        self.loads = []

    def __repr__(self):
        return("\nQ Load\n"
               "Name: %s\n"
               "Start Point: %s\n"
               "End Point: %s\n"
               "Start Value: %s\n"
               "End Value: %s\n"
               "Length: %s\n"
               "Total Load: %s\n" % (self.name,
                                     self.start_point,
                                     self.end_point,
                                     round(self.start_val, 2),
                                     round(self.end_val, 2),
                                     round(self.L, 2),
                                     round(self.total_load, 2)))

    def __str__(self):
        return("\nQ Load\n"
               "Name: %s\n"
               "Start Point: %s\n"
               "End Point: %s\n"
               "Start Value: %s\n"
               "End Value: %s\n"
               "Length: %s\n"
               "Total Load: %s\n" % (self.name,
                                     self.start_point,
                                     self.end_point,
                                     round(self.start_val, 2),
                                     round(self.end_val, 2),
                                     round(self.L, 2),
                                     round(self.total_load, 2)))

    def process_nodes(self, system):
        """
        This method distributes the distributed load to the nodes in the
        structural system according to the distribution method in the
        method distribute.

        If the load is shorter than a beam or ends halfway a beam in a point
        where no node is present, a new node is created in the end point of the
        load.
        """
        s = system
        p_s = self.start_point
        p_e = self.end_point

        if self.beam is None:

            l_el = len(s.elements)

            for e in range(l_el):

                element = s.elements[e]

                n_1 = element.node_1.point
                n_2 = element.node_2.point

                if n_1.on_line(p_s, p_e):

                    if n_2.on_line(p_s, p_e):

                        self.element_list.append(element)

                        q_0 = (self.s_load + self.d_load *
                               n_1.on_line(p_s, p_e, length=True)[1])
                        q_1 = (self.s_load + self.d_load *
                               n_2.on_line(p_s, p_e, length=True)[1])

                        T = element.trans_matrix

                        F = self.distribute(q_0, q_1, element.L, T)

                        self.loads.append(s.add_p_load(element.node_1,
                                                       name=self.name,
                                                       Fx=F[0],
                                                       Fy=F[1],
                                                       Fz=F[2],
                                                       Tx=F[3],
                                                       Ty=F[4],
                                                       Tz=F[5]))

                        self.loads.append(s.add_p_load(element.node_2,
                                                       name=self.name,
                                                       Fx=F[6],
                                                       Fy=F[7],
                                                       Fz=F[8],
                                                       Tx=F[9],
                                                       Ty=F[10],
                                                       Tz=F[11]))

                    else:

                        if np.all(p_e.on_line(n_1, n_2, length=True)):

                            ne = element.split_by_point(p_e, s)

                            self.element_list.append(element)

                            q_0 = (self.s_load + self.d_load *
                                   n_1.on_line(p_s, p_e, length=True)[1])
                            q_1 = self.e_load

                            T = element.trans_matrix

                            F = self.distribute(q_0, q_1, element.L, T)

                            self.loads.append(s.add_p_load(element.node_1,
                                                           name=self.name,
                                                           Fx=F[0],
                                                           Fy=F[1],
                                                           Fz=F[2],
                                                           Tx=F[3],
                                                           Ty=F[4],
                                                           Tz=F[5]))

                            self.loads.append(s.add_p_load(element.node_2,
                                                           name=self.name,
                                                           Fx=F[6],
                                                           Fy=F[7],
                                                           Fz=F[8],
                                                           Tx=F[9],
                                                           Ty=F[10],
                                                           Tz=F[11]))

                        elif np.all(p_s.on_line(n_1, n_2, length=True)):

                            ne = element.split_by_point(p_s, s)

                            self.element_list.append(element)

                            q_0 = (self.s_load + self.d_load *
                                   n_1.on_line(p_s, p_e, length=True)[1])
                            q_1 = self.s_load

                            T = element.trans_matrix

                            F = self.distribute(q_0, q_1, element.L, T)

                            self.loads.append(s.add_p_load(element.node_1,
                                                           name=self.name,
                                                           Fx=F[0],
                                                           Fy=F[1],
                                                           Fz=F[2],
                                                           Tx=F[3],
                                                           Ty=F[4],
                                                           Tz=F[5]))

                            self.loads.append(s.add_p_load(element.node_2,
                                                           name=self.name,
                                                           Fx=F[6],
                                                           Fy=F[7],
                                                           Fz=F[8],
                                                           Tx=F[9],
                                                           Ty=F[10],
                                                           Tz=F[11]))

                else:

                    if n_2.on_line(p_s, p_e):

                        if np.all(p_s.on_line(n_1, n_2, length=True)):

                            ne = element.split_by_point(p_s, s)
                            self.element_list.append(ne)

                            q_0 = self.s_load
                            q_1 = (self.s_load + self.d_load *
                                   n_2.on_line(p_s, p_e, length=True)[1])

                            T = ne.trans_matrix

                            F = self.distribute(q_0, q_1, ne.L, T)
                            self.loads.append(s.add_p_load(ne.node_1,
                                                           name=self.name,
                                                           Fx=F[0],
                                                           Fy=F[1],
                                                           Fz=F[2],
                                                           Tx=F[3],
                                                           Ty=F[4],
                                                           Tz=F[5]))

                            self.loads.append(s.add_p_load(ne.node_2,
                                                           name=self.name,
                                                           Fx=F[6],
                                                           Fy=F[7],
                                                           Fz=F[8],
                                                           Tx=F[9],
                                                           Ty=F[10],
                                                           Tz=F[11]))

                        elif np.all(p_e.on_line(n_1, n_2, length=True)):

                            ne = element.split_by_point(p_e, s)
                            self.element_list.append(ne)

                            q_0 = self.e_load
                            q_1 = (self.s_load + self.d_load *
                                   n_2.on_line(p_s, p_e, length=True)[1])

                            T = ne.trans_matrix

                            F = self.distribute(q_0, q_1, ne.L, T)

                            self.loads.append(s.add_p_load(ne.node_1,
                                                           name=self.name,
                                                           Fx=F[0],
                                                           Fy=F[1],
                                                           Fz=F[2],
                                                           Tx=F[3],
                                                           Ty=F[4],
                                                           Tz=F[5]))

                            self.loads.append(s.add_p_load(ne.node_2,
                                                           name=self.name,
                                                           Fx=F[6],
                                                           Fy=F[7],
                                                           Fz=F[8],
                                                           Tx=F[9],
                                                           Ty=F[10],
                                                           Tz=F[11]))

        else:

            l_el = len(self.beam.elements)

            for e in range(l_el):

                element = self.beam.elements[e]

                n_1 = element.node_1.point
                n_2 = element.node_2.point

                if n_1.on_line(p_s, p_e):

                    if n_2.on_line(p_s, p_e):

                        self.element_list.append(element)

                        q_0 = (self.s_load + self.d_load *
                               n_1.on_line(p_s, p_e, length=True)[1])
                        q_1 = (self.s_load + self.d_load *
                               n_2.on_line(p_s, p_e, length=True)[1])

                        T = element.trans_matrix

                        F = self.distribute(q_0, q_1, element.L, T)

                        self.loads.append(s.add_p_load(element.node_1,
                                                       name=self.name,
                                                       Fx=F[0],
                                                       Fy=F[1],
                                                       Fz=F[2],
                                                       Tx=F[3],
                                                       Ty=F[4],
                                                       Tz=F[5]))

                        self.loads.append(s.add_p_load(element.node_2,
                                                       name=self.name,
                                                       Fx=F[6],
                                                       Fy=F[7],
                                                       Fz=F[8],
                                                       Tx=F[9],
                                                       Ty=F[10],
                                                       Tz=F[11]))

                    else:

                        if np.all(p_e.on_line(n_1, n_2, length=True)):

                            ne = element.split_by_point(p_e, s)

                            self.element_list.append(element)

                            q_0 = (self.s_load + self.d_load *
                                   n_1.on_line(p_s, p_e, length=True)[1])
                            q_1 = self.e_load

                            T = element.trans_matrix

                            F = self.distribute(q_0, q_1, element.L, T)

                            self.loads.append(s.add_p_load(element.node_1,
                                                           name=self.name,
                                                           Fx=F[0],
                                                           Fy=F[1],
                                                           Fz=F[2],
                                                           Tx=F[3],
                                                           Ty=F[4],
                                                           Tz=F[5]))

                            self.loads.append(s.add_p_load(element.node_2,
                                                           name=self.name,
                                                           Fx=F[6],
                                                           Fy=F[7],
                                                           Fz=F[8],
                                                           Tx=F[9],
                                                           Ty=F[10],
                                                           Tz=F[11]))

                        elif np.all(p_s.on_line(n_1, n_2, length=True)):

                            ne = element.split_by_point(p_s, s)

                            self.element_list.append(element)

                            q_0 = (self.s_load + self.d_load *
                                   n_1.on_line(p_s, p_e, length=True)[1])
                            q_1 = self.s_load

                            T = element.trans_matrix

                            F = self.distribute(q_0, q_1, element.L, T)

                            self.loads.append(s.add_p_load(element.node_1,
                                                           name=self.name,
                                                           Fx=F[0],
                                                           Fy=F[1],
                                                           Fz=F[2],
                                                           Tx=F[3],
                                                           Ty=F[4],
                                                           Tz=F[5]))

                            self.loads.append(s.add_p_load(element.node_2,
                                                           name=self.name,
                                                           Fx=F[6],
                                                           Fy=F[7],
                                                           Fz=F[8],
                                                           Tx=F[9],
                                                           Ty=F[10],
                                                           Tz=F[11]))

                else:

                    if n_2.on_line(p_s, p_e):

                        if np.all(p_s.on_line(n_1, n_2, length=True)):

                            ne = element.split_by_point(p_s, s)
                            self.element_list.append(ne)

                            q_0 = self.s_load
                            q_1 = (self.s_load + self.d_load *
                                   n_2.on_line(p_s, p_e, length=True)[1])

                            T = ne.trans_matrix

                            F = self.distribute(q_0, q_1, ne.L, T)

                            self.loads.append(s.add_p_load(ne.node_1,
                                                           name=self.name,
                                                           Fx=F[0],
                                                           Fy=F[1],
                                                           Fz=F[2],
                                                           Tx=F[3],
                                                           Ty=F[4],
                                                           Tz=F[5]))

                            self.loads.append(s.add_p_load(ne.node_2,
                                                           name=self.name,
                                                           Fx=F[6],
                                                           Fy=F[7],
                                                           Fz=F[8],
                                                           Tx=F[9],
                                                           Ty=F[10],
                                                           Tz=F[11]))

                        elif np.all(p_e.on_line(n_1, n_2, length=True)):

                            ne = element.split_by_point(p_e, s)
                            self.element_list.append(ne)

                            q_0 = self.e_load
                            q_1 = (self.s_load + self.d_load *
                                   n_2.on_line(p_s, p_e, length=True)[1])

                            T = ne.trans_matrix

                            F = self.distribute(q_0, q_1, ne.L, T)

                            self.loads.append(s.add_p_load(ne.node_1,
                                                           name=self.name,
                                                           Fx=F[0],
                                                           Fy=F[1],
                                                           Fz=F[2],
                                                           Tx=F[3],
                                                           Ty=F[4],
                                                           Tz=F[5]))

                            self.loads.append(s.add_p_load(ne.node_2,
                                                           name=self.name,
                                                           Fx=F[6],
                                                           Fy=F[7],
                                                           Fz=F[8],
                                                           Tx=F[9],
                                                           Ty=F[10],
                                                           Tz=F[11]))

    def distribute(self, q_0_arr, q_1_arr, L, T):
        """
        This method converts the loads in the local axis system of the node
        they are applied to.
        Then it calculates the moments and shear forces to be applied to the
        nodes in their local axis system.
        It returns the loads in the global axis system to be applied to the
        nodes.
        """

        F_vect = np.zeros(12)

        q_0_arr = np.dot(T[:q_0_arr.shape[0], :q_0_arr.shape[0]], q_0_arr)
        q_1_arr = np.dot(T[:q_1_arr.shape[0], :q_1_arr.shape[0]], q_1_arr)

        for q_0, q_1, i in zip(q_0_arr, q_1_arr, range(q_0_arr.shape[0])):

            V1 = (0.35 * q_0 + 0.15 * q_1) * L
            M1 = ((1.5 * q_0 + q_1) / 30) * L**2
            V2 = (0.35 * q_1 + 0.15 * q_0) * L
            M2 = -((1.5 * q_1 + q_0) / 30) * L**2

            if np.isclose(V1+V2, q_0*L + (q_1-q_0)*L/2):
                if i < 3:
                    F_vect[i] = V1
                    F_vect[i+6] = V2
                    if i != 0:
                        F_vect[6-i] = M1
                        F_vect[12-i] = M2
            else:
                raise Exception("Distributed load cannot be distributed " +
                                "properly at " + self)

        F = np.dot(T.T, F_vect)

        return F


class SquareLoad():
    """
    This class contains a rectangular load, which distributes the constant
    load value supplied over the bordering beams.
    It has to be bordered by the beams on which the load is to be
    distributed in its current form.
    """

    def __init__(self, value, start_point, end_point,
                 direction=[0, 0, -1], beams=None, name='sq load'):
        """
        Initializes the SquareLoad class.

        value: Constant value of the load (N/mm²)
        start_point: Point class containing the start location of the rectangle
        end_point: Point class containing the end location of the rectangle
        direction: Vector containg the direction of the load
        beams: Beam class objects the load is constrained to
        name: The name of the applied load
        """

        self.name = name
        self.q = value
        self.start_point = start_point
        self.end_point = end_point
        self.direction = direction
        self.beams = beams

        self.sys = None

        self.third_point, self.fourth_point, self.sides = self.complete_points(
                                              self.start_point, self.end_point)

        self.A = self.sides[0]*self.sides[1]
        self.total_load = np.abs(self.A * self.q)

        self.sq_loads = None

    def __repr__(self):
        return("\nRectangular Load\n"
               "Name: %s\n"
               "Start Point: %s\n"
               "End Point: %s\n"
               "Value: %s\n"
               "Area: %s\n"
               "Total Load: %s\n" % (self.name,
                                     self.start_point,
                                     self.end_point,
                                     round(self.q, 2),
                                     round(self.A, 2),
                                     round(self.total_load, 2)))

    def __str__(self):
        return("\nRectangular Load\n"
               "Name: %s\n"
               "Start Point: %s\n"
               "End Point: %s\n"
               "Value: %s\n"
               "Area: %s\n"
               "Total Load: %s\n" % (self.name,
                                     self.start_point,
                                     self.end_point,
                                     round(self.q, 2),
                                     round(self.A, 2),
                                     round(self.total_load, 2)))

    def complete_points(self, p_s, p_e):
        """
        Calculates the positions of points 3 and 4 of the rectangle.
        """

        ps = p_e.v - p_s.v

        sides = ps[np.abs(ps) > 0]

        try:
            ps_1 = np.where(ps == sides[0], sides[0], 0)
            ps_2 = np.where(ps == sides[1], sides[1], 0)
        except IndexError:
            print()
            print(p_s)
            print(p_e)
            print(ps)

        third_point = Point(*(p_s.v + ps_2))
        fourth_point = Point(*(p_s.v + ps_1))

        return third_point, fourth_point, sides

    def process_beams(self, system=None):
        """
        Checks if there are open sides to the square and returns the points
        comprising these open sides.
        """

        if system is None:
            if self.sys is None:
                raise ValueError('Please give the square load a system to ' +
                                 'work in')
        else:
            self.sys = system

        p_s = self.start_point
        p_t = self.third_point
        p_e = self.end_point
        p_f = self.fourth_point

        ps_pt = False
        ps_pf = False
        pt_pe = False
        pf_pe = False

        self.beams = []

        for beam in self.sys.beams:

            if self.check_beam(beam, p_s, p_t):
                ps_pt = True
                self.beams.append(beam)
            elif self.check_beam(beam, p_s, p_f):
                ps_pf = True
                self.beams.append(beam)
            elif self.check_beam(beam, p_t, p_e):
                pt_pe = True
                self.beams.append(beam)
            elif self.check_beam(beam, p_f, p_e):
                pf_pe = True
                self.beams.append(beam)

        s_o = []

        if ps_pt is False:
            if p_s not in s_o:
                s_o.append(p_s)
            if p_f not in s_o:
                s_o.append(p_t)

        if ps_pf is False:
            if p_s not in s_o:
                s_o.append(p_s)
            if p_f not in s_o:
                s_o.append(p_f)

        if pt_pe is False:
            if p_t not in s_o:
                s_o.append(p_t)
            if p_e not in s_o:
                s_o.append(p_e)

        if pf_pe is False:
            if p_f not in s_o:
                s_o.append(p_f)
            if p_e not in s_o:
                s_o.append(p_e)

        if len(s_o) == 4:
            s_o.insert(2, s_o[0].in_between(s_o[1]))
            s_o.append(s_o[3].in_between(s_o[4]))

        return s_o

    def process_beams2(self, system=None):
        """
        Looks to collect all beams within the square load square - subsequently
        looking to form subrectangles to distribute load to.
        """

        if system is None:
            if self.sys is None:
                raise ValueError('Please give the square load a system to ' +
                                 'work in')
        else:
            self.sys = system

        self.sq_loads = []

        p_s = self.start_point
        p_t = self.third_point
        p_e = self.end_point
        p_f = self.fourth_point

        if self.beams is None:
            self.beams = []

        p_list = []

        for beam in self.sys.beams:

            check = self.check_beam_rect(beam, p_s, p_f, p_e, p_t)

            if check:

                self.beams.append(beam)

                p_list.append([check[1], check[2]])

        self.rectangle = Rectangle(p_s, p_e, p_f, p_t)

        p_list = [p for k in p_list for p in k]

        p_list2 = []
        for p in p_list:
            if p not in p_list2:
                p_list2.append(p)

        p_list = p_list2

        self.rectangle.find_children(p_list)

        if self.rectangle.children:
            r_list = self.rectangle.children
            if self.rectangle.c_children is not None:
                r_list += self.rectangle.c_children
            r_list = [r for r in r_list if r.smallest]
        else:
            r_list = [self.rectangle]

        for rect in r_list:

            bm = [b for b in self.beams if
                  self.check_beam(b, rect.A, rect.B)]
            bm += [b for b in self.beams if
                   self.check_beam(b, rect.B, rect.C)]
            bm += [b for b in self.beams if
                   self.check_beam(b, rect.C, rect.D)]
            bm += [b for b in self.beams if
                   self.check_beam(b, rect.D, rect.A)]

            self.sq_loads.append(SquareLoad(self.q,
                                            rect.A,
                                            rect.C,
                                            direction=self.direction,
                                            beams=bm,
                                            name='auto sq load'))

            s_o = self.sq_loads[-1].process_beams(self.sys)
            self.sq_loads[-1].distribute_load(s_o)

    def check_beam(self, beam, p_1, p_2):
        """
        Checks if beam is contained by the line between p_1 and p_2.
        If it does not, checks if p_1 or p_2 are contained by the line
        between the beams end nodes.
        """
        if beam.sn_1.point.on_line(p_1, p_2):
            if beam.sn_2.point.on_line(p_1, p_2):
                return True
            else:
                if p_1.on_line(beam.sn_1.point, beam.sn_2.point):
                    if p_2.on_line(beam.sn_1.point, beam.sn_2.point):
                        return True
                    else:
                        return False
                else:
                    return False
        else:
            if p_1.on_line(beam.sn_1.point, beam.sn_2.point):
                if p_2.on_line(beam.sn_1.point, beam.sn_2.point):
                    return True
                else:
                    return False
            else:
                return False

    def check_beam_rect(self, beam, p_1, p_2, p_3, p_4):
        """
        Checks if a (portion of a) beam is inside of the given rectangle.
        If True, returns the points comprising the portion of the beam
        inside the rectangle.
        """
        pb_1 = beam.sn_1.point
        pb_2 = beam.sn_2.point
        if pb_1.in_rect(p_1, p_2, p_3, p_4):
            if pb_2.in_rect(p_1, p_2, p_3, p_4):
                return True, pb_1, pb_2
            else:
                if pb_1.intersect(pb_1, pb_2, p_1, p_2):
                    drop, pr_2 = pb_1.intersect(pb_1, pb_2, p_1, p_2)
                    return True, pb_1, Point(*pr_2)
                elif pb_1.intersect(pb_1, pb_2, p_2, p_3):
                    drop, pr_2 = pb_1.intersect(pb_1, pb_2, p_2, p_3)
                    return True, pb_1, Point(*pr_2)
                elif pb_1.intersect(pb_1, pb_2, p_3, p_4):
                    drop, pr_2 = pb_1.intersect(pb_1, pb_2, p_3, p_4)
                    return True, pb_1, Point(*pr_2)
                elif pb_1.intersect(pb_1, pb_2, p_4, p_1):
                    drop, pr_2 = pb_1.intersect(pb_1, pb_2, p_4, p_1)
                    return True, pb_1, Point(*pr_2)
                else:
                    return False
        elif pb_2.in_rect(p_1, p_2, p_3, p_4):
            if pb_2.intersect(pb_1, pb_2, p_1, p_2):
                drop, pr_2 = pb_1.intersect(pb_1, pb_2, p_1, p_2)
                return True, pb_2, Point(*pr_2)
            elif pb_2.intersect(pb_1, pb_2, p_2, p_3):
                drop, pr_2 = pb_1.intersect(pb_1, pb_2, p_2, p_3)
                return True, pb_2, Point(*pr_2)
            elif pb_2.intersect(pb_1, pb_2, p_3, p_4):
                drop, pr_2 = pb_1.intersect(pb_1, pb_2, p_3, p_4)
                return True, pb_2, Point(*pr_2)
            elif pb_2.intersect(pb_1, pb_2, p_4, p_1):
                drop, pr_2 = pb_1.intersect(pb_1, pb_2, p_4, p_1)
                return True, pb_2, Point(*pr_2)
            else:
                return False
        else:
            return False

    def distribute_load(self, s_o=None):
        """
        Distributes the loads to the beams on four, three or two sides of the
        square.

        https://engineering.stackexchange.com/questions/5638/
        how-to-calculate-the-forces-on-members-of-a-frame-due-to-a-
        distributed-pressure
        """

        if self.sys is None:
            raise AttributeError("System not defined for square load " +
                                 "distribution")

        if s_o:

            if len(s_o) == 2:

                self.one_side_open(s_o)

            elif len(s_o) == 3:

                self.two_con_sides_open(s_o)

            elif len(s_o) == 4:
                raise ValueError("Can't determine open sides from four " +
                                 "points")
            elif len(s_o) == 6:

                self.two_opposing_sides_open(s_o)

        else:

            self.all_closed()

    def one_side_open(self, s_o):
        """
        This method distributes the load to three closed sides of the square.
        """

        q = self.q
        d = self.direction
        p_s = self.start_point
        p_e = self.end_point

        load_c = []

        s_o_0 = s_o[0]
        s_o_1 = s_o[1]

        s_2 = np.max(np.abs(s_o_1.v-s_o_0.v))
        p_2 = [k for k in [p_s, p_e] if k not in s_o][0]

        s_o_2 = [p for p in [s_o_0.v, s_o_1.v] if
                 np.count_nonzero(p - p_2.v) < 2][0]

        s_1 = np.abs(s_o_2 - p_2.v)[np.abs(s_o_2 - p_2.v) > 0][0]

        if s_2 / 2 > s_1:

            d_2 = np.where(np.abs(s_o_1.v-s_o_0.v) == s_2)[0]
            d_ = np.delete(np.arange(3), d_2)

            if p_s.v[d_2] < p_e.v[d_2]:

                l_2 = np.array([p_s.v[d_2],
                                p_s.v[d_2] + s_1,
                                p_e.v[d_2] - s_1,
                                p_e.v[d_2]])

                f_2 = np.array([0, q * s_1, q * s_1, 0])

            else:

                l_2 = np.array([p_e.v[d_2],
                                p_e.v[d_2] + s_1,
                                p_s.v[d_2] - s_1,
                                p_s.v[d_2]])

                f_2 = np.array([0, q * s_1, q * s_1, 0])

        else:

            d_2 = np.where(np.abs(s_o_1.v-s_o_0.v) == s_2)[0]
            d_ = np.delete(np.arange(3), d_2)

            if p_s.v[d_2] < p_e.v[d_2]:

                l_2 = np.linspace(p_s.v[d_2], p_e.v[d_2], 3)

            else:

                l_2 = np.linspace(p_e.v[d_2], p_s.v[d_2], 3)

            f_2 = np.array([0, q * s_2 / 2, 0])

        for c in range(1, len(l_2)):

            xyz = np.zeros(3)
            xyz[d_] = p_2.v[d_]
            xyz[d_2] = l_2[c-1]
            p_l_1 = Point(*xyz)
            xyz[d_2] = l_2[c]
            p_l_2 = Point(*xyz)

            if self.beams is not None:

                try:

                    bm = [b for b in self.beams if
                          self.check_beam(b, p_l_1, p_l_2)][0]

                    load_c.append(self.sys.add_q_load(f_2[c-1],
                                                      p_l_1,
                                                      p_l_2,
                                                      end_val=f_2[c],
                                                      beam=bm,
                                                      direction=d,
                                                      name=self.name
                                                      ).total_load)
                except IndexError:

                    load_c.append(self.sys.add_q_load(f_2[c-1],
                                                      p_l_1,
                                                      p_l_2,
                                                      end_val=f_2[c],
                                                      direction=d,
                                                      name=self.name
                                                      ).total_load)
            else:

                load_c.append(self.sys.add_q_load(f_2[c-1],
                                                  p_l_1,
                                                  p_l_2,
                                                  end_val=f_2[c],
                                                  direction=d,
                                                  name=self.name
                                                  ).total_load)

        if s_1 == s_2:

            d_1 = np.where(np.abs(s_o_2 - p_2.v) == s_1)[0]
            d_ = np.delete(np.arange(3), d_1)
            l_1 = np.array([s_o_2[d_1], s_o_2[d_1]+l_2[1], p_2.v[d_1]])

            if s_2 > s_1:

                f_1 = np.array([0, q * s_1 / 2, q * s_1 / 2])

            else:

                f_1 = np.array([0, q * s_2 / 2, q * s_2 / 2])

            for c in range(1, len(l_1)):

                for point in [s_o_0, s_o_1]:

                    xyz = np.zeros(3)
                    xyz[d_] = point.v[d_]
                    xyz[d_1] = l_1[c-1]
                    p_l_1 = Point(*xyz)
                    xyz[d_1] = l_1[c]
                    p_l_2 = Point(*xyz)

                    if self.beams is not None:

                        try:

                            bm = [b for b in self.beams if
                                  self.check_beam(b, p_l_1, p_l_2)][0]

                            load_c.append(self.sys.add_q_load(f_1[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_1[c],
                                                              direction=d,
                                                              beam=bm,
                                                              name=self.name
                                                              ).total_load)

                        except IndexError:

                            load_c.append(self.sys.add_q_load(f_1[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_1[c],
                                                              direction=d,
                                                              beam=bm,
                                                              name=self.name
                                                              ).total_load)
                    else:

                        load_c.append(self.sys.add_q_load(f_1[c-1],
                                                          p_l_1,
                                                          p_l_2,
                                                          end_val=f_1[c],
                                                          direction=d,
                                                          name=self.name
                                                          ).total_load)
        else:

            d_1 = np.where(np.abs(s_o_2 - p_2.v) == s_1)[0]
            d_ = np.delete(np.arange(3), d_1)

            if s_2 / 2 > s_1:

                f_1 = np.array([0, q * s_1])

                if s_o_2[d_1] > p_2.v[d_1]:

                    l_1 = np.array([p_2.v[d_1],
                                    s_o_2[d_1]])

                else:

                    l_1 = np.array([s_o_2[d_1],
                                    p_2.v[d_1]])

            else:

                if s_o_2[d_1] > p_2.v[d_1]:

                    l_1 = np.array([p_2.v[d_1],
                                    p_2.v[d_1]+(l_2[1]-l_2[0]),
                                    s_o_2[d_1]])

                else:

                    l_1 = np.array([s_o_2[d_1],
                                    s_o_2[d_1]+(l_2[1]-l_2[0]),
                                    p_2.v[d_1]])

                f_1 = np.array([0, q * s_2 / 2, q * s_2 / 2, q * s_2 / 2])

            for c in range(1, len(l_1)):

                for point in [s_o_0, s_o_1]:

                    xyz = np.zeros(3)
                    xyz[d_] = point.v[d_]
                    xyz[d_1] = l_1[c-1]
                    p_l_1 = Point(*xyz)
                    xyz[d_1] = l_1[c]
                    p_l_2 = Point(*xyz)

                    if self.beams is not None:

                        try:

                            bm = [b for b in self.beams if
                                  self.check_beam(b, p_l_1, p_l_2)][0]

                            load_c.append(self.sys.add_q_load(f_1[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_1[c],
                                                              direction=d,
                                                              beam=bm,
                                                              name=self.name
                                                              ).total_load)

                        except IndexError:

                            load_c.append(self.sys.add_q_load(f_1[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_1[c],
                                                              direction=d,
                                                              name=self.name
                                                              ).total_load)

                    else:

                        load_c.append(self.sys.add_q_load(f_1[c-1],
                                                          p_l_1,
                                                          p_l_2,
                                                          end_val=f_1[c],
                                                          direction=d,
                                                          name=self.name
                                                          ).total_load)

            load_check = np.sum(load_c) - self.total_load

            if not np.isclose(load_check, 0):
                raise Exception("Total square load not equal to imposed load")

    def two_con_sides_open(self, s_o):
        """
        This method distributes the load on the two consecutive closed sides
        of the square.
        """

        q = self.q
        d = self.direction
        p_s = self.start_point
        p_t = self.third_point
        p_e = self.end_point
        p_f = self.fourth_point

        p_list = [p_s, p_t, p_e, p_f]

        p_cc = [p for p in p_list if p not in s_o][0]

        p_cn = [p for p in p_list if
                np.equal(p, p_cc)[np.equal(p, p_cc) is True].shape[0] == 2]

        s_1 = np.max(np.abs(p_cn[0] - p_cc))
        s_2 = np.max(np.abs(p_cn[1] - p_cc))

        q_1 = [0, q*(s_2/s_1)]
        q_2 = [0, q*(s_1/s_2)]

        load_check = 0

        if self.beams is not None:

            try:

                bm_1 = [b for b in self.beams if
                        self.check_beam(b, p_cc, p_cn[0])][0]

                load_check += self.sys.add_q_load(q_1[0],
                                                  p_cc,
                                                  p_cn[0],
                                                  end_val=q_1[1],
                                                  direction=d,
                                                  beam=bm_1,
                                                  name=self.name
                                                  ).total_load

            except IndexError:

                load_check += self.sys.add_q_load(q_1[0],
                                                  p_cc,
                                                  p_cn[0],
                                                  end_val=q_1[1],
                                                  direction=d,
                                                  name=self.name
                                                  ).total_load

            try:

                bm_2 = [b for b in self.beams if
                        self.check_beam(b, p_cc, p_cn[1])][0]

                load_check += self.sys.add_q_load(q_2[0],
                                                  p_cc,
                                                  p_cn[1],
                                                  end_val=q_2[1],
                                                  direction=d,
                                                  beam=bm_2,
                                                  name=self.name
                                                  ).total_load

            except IndexError:

                load_check += self.sys.add_q_load(q_2[0],
                                                  p_cc,
                                                  p_cn[1],
                                                  end_val=q_2[1],
                                                  direction=d,
                                                  name=self.name
                                                  ).total_load

        else:

            load_check += self.sys.add_q_load(q_1[0],
                                              p_cc,
                                              p_cn[0],
                                              end_val=q_1[1],
                                              direction=d,
                                              name=self.name
                                              ).total_load

            load_check += self.sys.add_q_load(q_2[0],
                                              p_cc,
                                              p_cn[1],
                                              end_val=q_2[1],
                                              direction=d,
                                              name=self.name
                                              ).total_load

        load_check -= self.total_load

        if not np.isclose(load_check, 0):
            raise Exception("Total square load not equal to imposed load")

    def two_opposing_sides_open(self, s_o):
        """
        This method distributes the load on the two opposing closed sides of
        the square
        """

        q = self.q
        d = self.direction

        load_check = 0

        if s_o[2].on_line(s_o[0], s_o[1]):
            if s_o[5].on_line(s_o[3], s_o[4]):
                s_q = np.max(np.abs(s_o[3].v-s_o[0].v))
                q *= s_q / 2

                if self.beams is not None:

                    try:

                        bm_1 = [b for b in self.beams if
                                self.check_beam(b, s_o[0], s_o[1])][0]

                        load_check += self.sys.add_q_load(q,
                                                          s_o[0],
                                                          s_o[1],
                                                          direction=d,
                                                          beam=bm_1,
                                                          name=self.name
                                                          ).total_load
                    except IndexError:

                        load_check += self.sys.add_q_load(q,
                                                          s_o[0],
                                                          s_o[1],
                                                          direction=d,
                                                          name=self.name
                                                          ).total_load

                    try:

                        bm_2 = [b for b in self.beams if
                                self.check_beam(b, s_o[3], s_o[4])][0]

                        load_check += self.sys.add_q_load(q,
                                                          s_o[3],
                                                          s_o[4],
                                                          direction=d,
                                                          beam=bm_2,
                                                          name=self.name
                                                          ).total_load

                    except IndexError:

                        load_check += self.sys.add_q_load(q,
                                                          s_o[3],
                                                          s_o[4],
                                                          direction=d,
                                                          name=self.name
                                                          ).total_load
                else:

                    load_check += self.sys.add_q_load(q,
                                                      s_o[0],
                                                      s_o[1],
                                                      direction=d,
                                                      name=self.name
                                                      ).total_load

                    load_check += self.sys.add_q_load(q,
                                                      s_o[3],
                                                      s_o[4],
                                                      direction=d,
                                                      name=self.name
                                                      ).total_load

        load_check -= self.total_load

        if not np.isclose(load_check, 0):
            raise Exception("Total square load not equal to imposed load")

    def all_closed(self):
        """
        This method distributes the load to the four closed sides of the square
        """

        q = self.q
        d = self.direction
        p_s = self.start_point
        p_e = self.end_point

        s_1 = np.max(np.abs(p_e.v-p_s.v))
        s_2 = np.min(np.abs(p_e.v-p_s.v)[np.abs(p_e.v-p_s.v) > 0])

        if s_1 == s_2:

            d_1 = np.where(np.abs(p_e.v-p_s.v) == s_2)[0]
            d_2 = np.where(np.abs(p_e.v-p_s.v) == s_2)[0]
            d_1_ = np.delete(np.arange(3), d_1)
            d_2_ = np.delete(np.arange(3), d_2)
            l_1 = np.linspace(p_s.v[d_1], p_e.v[d_1], 3)
            l_2 = np.linspace(p_s.v[d_2], p_e.v[d_2], 3)
            f_12 = np.array([0, q * s_2 / 2, 0])

            load_check = 0

            for c in range(1, len(l_1)):

                for point in [p_s, p_e]:

                    xyz = np.zeros(3)
                    xyz[d_1_] = point.v[d_1]
                    xyz[d_1] = l_1[c-1]
                    p_l_1 = Point(*xyz)
                    xyz[d_1] = l_1[c]
                    p_l_2 = Point(*xyz)

                    if self.beams is not None:

                        try:

                            bm_1 = [b for b in self.beams if
                                    self.check_beam(b, p_l_1, p_l_2)][0]

                            load_check += self.sys.add_q_load(f_12[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_12[c],
                                                              direction=d,
                                                              beam=bm_1,
                                                              name=self.name
                                                              ).total_load
                        except IndexError:

                            load_check += self.sys.add_q_load(f_12[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_12[c],
                                                              direction=d,
                                                              name=self.name
                                                              ).total_load
                    else:

                        load_check += self.sys.add_q_load(f_12[c-1], p_l_1,
                                                          p_l_2,
                                                          end_val=f_12[c],
                                                          direction=d,
                                                          name=self.name
                                                          ).total_load

                    xyz = np.zeros(3)
                    xyz[d_2_] = point.v[d_2]
                    xyz[d_2] = l_2[c-1]
                    p_l_1 = Point(*xyz)
                    xyz[d_2] = l_2[c]
                    p_l_2 = Point(*xyz)

                    if self.beams is not None:

                        try:

                            bm_2 = [b for b in self.beams if
                                    self.check_beam(b, p_l_1, p_l_2)][0]

                            load_check += self.sys.add_q_load(f_12[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_12[c],
                                                              direction=d,
                                                              beam=bm_2,
                                                              name=self.name
                                                              ).total_load

                        except IndexError:

                            load_check += self.sys.add_q_load(f_12[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_12[c],
                                                              direction=d,
                                                              name=self.name
                                                              ).total_load

                    else:

                        load_check += self.sys.add_q_load(f_12[c-1],
                                                          p_l_1,
                                                          p_l_2,
                                                          end_val=f_12[c],
                                                          direction=d,
                                                          name=self.name
                                                          ).total_load

            load_check -= self.total_load

            if not np.isclose(load_check, 0):
                raise Exception("Total square load not equal to imposed load")

        else:

            d_2 = np.where(np.abs(p_e.v-p_s.v) == s_2)[0]
            d_ = np.delete(np.arange(3), d_2)
            if p_s.v[d_2][0] > p_e.v[d_2][0]:
                l_2 = np.linspace(p_e.v[d_2], p_s.v[d_2], 3)
            else:
                l_2 = np.linspace(p_s.v[d_2], p_e.v[d_2], 3)
            f_2 = np.array([0, q * s_2 / 2, 0])

            load_check = 0

            for c in range(1, len(l_2)):

                for point in [p_s, p_e]:

                    xyz = np.zeros(3)
                    xyz[d_] = point.v[d_]
                    xyz[d_2] = l_2[c-1]
                    p_l_1 = Point(*xyz)
                    xyz[d_2] = l_2[c]
                    p_l_2 = Point(*xyz)

                    if self.beams is not None:

                        try:

                            bm_2 = [b for b in self.beams if
                                    self.check_beam(b, p_l_1, p_l_2)][0]

                            load_check += self.sys.add_q_load(f_2[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_2[c],
                                                              direction=d,
                                                              beam=bm_2,
                                                              name=self.name
                                                              ).total_load

                        except IndexError:

                            load_check += self.sys.add_q_load(f_2[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_2[c],
                                                              direction=d,
                                                              name=self.name
                                                              ).total_load

                    else:

                        load_check += self.sys.add_q_load(f_2[c-1],
                                                          p_l_1,
                                                          p_l_2,
                                                          end_val=f_2[c],
                                                          direction=d,
                                                          name=self.name
                                                          ).total_load

            d_1 = np.where(np.abs(p_e.v-p_s.v) == s_1)[0]
            d_ = np.delete(np.arange(3), d_1)
            if p_s.v[d_1][0] > p_e.v[d_1][0]:
                l_1 = np.array([p_e.v[d_1], p_e.v[d_1]+(l_2[1]-l_2[0]),
                                p_s.v[d_1]-(l_2[1]-l_2[0]), p_s.v[d_1]])
            else:
                l_1 = np.array([p_s.v[d_1], p_s.v[d_1]+(l_2[1]-l_2[0]),
                                p_e.v[d_1]-(l_2[1]-l_2[0]), p_e.v[d_1]])
            f_1 = np.array([0, q * s_2 / 2, q * s_2 / 2, 0])

            for c in range(1, len(l_1)):

                for point in [p_s, p_e]:

                    xyz = np.zeros(3)
                    xyz[d_] = point.v[d_]
                    xyz[d_1] = l_1[c-1]
                    p_l_1 = Point(*xyz)
                    xyz[d_1] = l_1[c]
                    p_l_2 = Point(*xyz)

                    if self.beams is not None:

                        try:

                            bm_1 = [b for b in self.beams if
                                    self.check_beam(b, p_l_1, p_l_2)][0]

                            load_check += self.sys.add_q_load(f_1[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_1[c],
                                                              direction=d,
                                                              beam=bm_1,
                                                              name=self.name
                                                              ).total_load
                        except IndexError:

                            load_check += self.sys.add_q_load(f_1[c-1],
                                                              p_l_1,
                                                              p_l_2,
                                                              end_val=f_1[c],
                                                              direction=d,
                                                              name=self.name
                                                              ).total_load
                    else:

                        load_check += self.sys.add_q_load(f_1[c-1],
                                                          p_l_1,
                                                          p_l_2,
                                                          end_val=f_1[c],
                                                          direction=d,
                                                          name=self.name
                                                          ).total_load

            load_check -= self.total_load

            if not np.isclose(load_check, 0):
                raise Exception("Total square load not equal to imposed load")
