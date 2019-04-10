from edsFEM.points import Point
from edsFEM.nodes import Node, StructureNode
from edsFEM.beams import Beam
from edsFEM.elements import EBElement
from edsFEM.supports import Support
from edsFEM.loads import QLoad, SquareLoad, PointLoad
from edsFEM.assemblers import BeamAssembler
from edsFEM.solvers import SimpleSolver
from edsFEM.postprocessors import BeamPostProcessor


class System():
    """
    The System class is the first class to be initialized when performing
    structural analysis with edsFEM.
    It contains all methods to construct, support, load, assemble and solve
    structural models.
    """

    def __init__(self, dims=3):
        """
        Initializes the System class.

        It is not possible at this time to to assign the number of available
        dimensions for the structural system using the dims argument as this
        functionality is not yet active.
        """

        self.dims = dims
        self.structurenodes = None
        self.beams = None
        self.nodes = None
        self.elements = None
        self.supports = None
        self.displacements = None
        self.square_loads = None
        self.q_loads = None
        self.p_loads = None
        self.m_loads = None

        self.assembler = None

        self.k_matrix = None
        self.f_vector = None
        self.u_vector = None

        self.k_reduced = None
        self.f_reduced = None
        self.u_reduced = None

    def assemble(self, version='Beam', selfweight=None, t_test=False):
        """
        This method assembles the reduced stiffness matrix and displacement
        vector to be solved for the system.
        In its current iteration it uses the BeamAssembler class to assemble
        the stiffness matrix.
        It is possible to automatically assign self-weight to the beams by
        supplying a load factor for the selfweight argument (float > 0).
        """

        if version == 'Beam':
            self.assembler = BeamAssembler(self)

        self.assembler.assemble(selfweight=selfweight, t_test=t_test)

    def solve(self, version='Simple', t_test=False):
        """
        This method solves the matrix equation posed by the reduced stiffness
        matrix and displacement vector assembled in the assemble method.
        In its current iteration it calls upon the SimpleSolver class to
        perform the solving of the system.
        """

        if version == 'Simple':
            self.solver = SimpleSolver(self)

        self.solver.solve(t_test=t_test)

    def postprocess(self, version='Beam'):
        """
        This method postprocesses the results of solving the system and assigns
        and calculates displacements, support reactions, member internal loads
        and stresses.
        In its current iteration it uses the BeamPostProcessor class to perform
        these tasks.
        """

        if version == 'Beam':
            self.pprocessor = BeamPostProcessor(self)

        self.pprocessor.assign_displacements()
        self.pprocessor.assign_forces()

    def add_beam(self, point_1, point_2, profile=None, E=None, A=None,
                 I_y=None, I_z=None, G=None, I_t=None, Wy_el=None,
                 Wz_el=None, global_cs=None, a_l=None, sys_eq=False,
                 Wy_pl=None, Wz_pl=None, release=None, segments=12):
        """
        This method adds a beam between the two points specified with either
        the properties of the specified profile class or the properties given
        in the arguments.

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
        """

        sn_1 = self.add_structurenode(point_1)
        sn_2 = self.add_structurenode(point_2)
        if profile is None:
            beam = Beam(sn_1, sn_2, self, E=E, A=A, I_y=I_y, I_z=I_z, G=G,
                        I_t=I_t, Wy_el=Wy_el, Wz_el=Wz_el,
                        global_cs=global_cs, a_l=a_l, sys_eq=sys_eq,
                        Wy_pl=Wy_pl, Wz_pl=Wz_pl, release=release,
                        segments=segments, profile=profile)
            return beam
        else:
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
            beam = Beam(sn_1, sn_2, self, E=E, A=A, I_y=I_y, I_z=I_z, G=G,
                        I_t=I_t, Wy_el=Wy_el, Wz_el=Wz_el,
                        global_cs=global_cs, a_l=a_l, sys_eq=sys_eq,
                        Wy_pl=Wy_pl, Wz_pl=Wz_pl, release=release,
                        segments=segments, profile=profile)
            return beam

    def add_structurenode(self, point):
        """
        This method adds a structurenode (and thus a normal node, too, if it
        does not yet exist) at the Point class specified in the point argument.
        The structurenode indicates the ends of a beam or column.
        """

        if self.structurenodes is not None:
            for ss in self.structurenodes:
                if ss.point == point:
                    return ss

            if self.nodes is not None:
                for node in self.nodes:
                    if node.point == point:
                        ss = StructureNode(len(self.structurenodes)+1,
                                           node=node)
                        self.structurenodes.append(ss)
                        return ss

                if self.beams is not None:
                    for beam in self.beams:
                        if beam.check_point(point):
                            node = beam.add_node(point)
                            ss = StructureNode(len(self.structurenodes)+1,
                                               node=node)
                            self.structurenodes.append(ss)
                            return ss

                self.add_node(point)
                ss = StructureNode(1, node=self.nodes[-1])
                self.structurenodes.append(ss)
                return ss

            else:

                if self.beams is not None:
                    for beam in self.beams:
                        if beam.check_point(point):
                            node = beam.add_node(point)
                            ss = StructureNode(len(self.structurenodes)+1,
                                               node=node)
                            self.structurenodes.append(ss)
                            return ss

                self.add_node(point)
                ss = StructureNode(len(self.structurenodes)+1,
                                   node=self.nodes[-1])
                self.structurenodes.append(ss)
                return ss

        else:

            if self.nodes is not None:
                for node in self.nodes:
                    if node.point == point:
                        ss = StructureNode(1, node=node)
                        self.structurenodes = [ss]
                        return ss

            else:

                self.add_node(point)
                ss = StructureNode(1, node=self.nodes[-1])
                self.structurenodes = [ss]
                return ss

    def add_node(self, position, Fx=0, Fy=0, Fz=0, Tx=0, Ty=0, Tz=0,
                 ux=0, uy=0, uz=0, phi_x=0, phi_y=0, phi_z=0,
                 support=None, elements=None):
        """
        This method adds a Node at the specified position.

        position: may be a list of coordinates; [x, y, z]
                  or a tuple of coordinates; (x, y, z)
                  or a Point class
        Through the Fx - Tz Arguments it is possible to add forces and moments
        in/around the principle global axes.
        Through the ux - phi_z arguments it is possible to add displacements
        and rotations in/around the principle global axes.
        In the support argument a Support class may be specified for the node.
        In the elements argument a list of elements attached to the node may
        be specified.
        """

        if type(position) is list:
            point = Point(position[0], position[1], position[2])
        elif type(position) is tuple:
            point = Point(position[0], position[1], position[2])
        elif type(position) is Point:
            point = position

        if self.nodes is not None:
            ID = len(self.nodes)+1
            node = Node(ID, point, Fx=Fx, Fy=Fy, Fz=Fz, Tx=Tx, Ty=Ty, Tz=Tz,
                        ux=ux, uy=uy, uz=uz, phi_x=phi_x, phi_y=phi_y,
                        phi_z=phi_z, support=support, elements=elements)
            self.nodes.append(node)
        else:
            ID = 1
            node = Node(ID, point, Fx=Fx, Fy=Fy, Fz=Fz, Tx=Tx, Ty=Ty, Tz=Tz,
                        ux=ux, uy=uy, uz=uz, phi_x=phi_x, phi_y=phi_y,
                        phi_z=phi_z, support=support, elements=elements)
            self.nodes = [node]
        return node

    def add_element(self, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                    node_1, node_2, global_cs=None, a_l=None, sys_eq=False,
                    Wy_pl=None, Wz_pl=None, release=None, beam=None):
        """
        This method adds and element between the specified nodes.

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
        beam: Beam class object the element is a part of
        """

        if self.elements is not None:

            ID = len(self.elements)+1

            element = EBElement(ID, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                                node_1, node_2, global_cs=global_cs,
                                a_l=a_l, sys_eq=sys_eq, Wy_pl=Wy_pl,
                                Wz_pl=Wz_pl, release=release, beam=beam)

            self.elements.append(element)

        else:

            ID = 1
            element = EBElement(ID, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el,
                                node_1, node_2, global_cs=global_cs,
                                a_l=a_l, sys_eq=sys_eq, Wy_pl=Wy_pl,
                                Wz_pl=Wz_pl, release=release, beam=beam)

            self.elements = [element]

        return element

    def add_support(self, node, t_x=False, t_y=False, t_z=False,
                    r_x=False, r_y=False, r_z=False):
        """
        This method adds a Support to the specified node.

        node: the node the support is attached to
        t_x: 0 - fixed
             any number != 0 - support stiffness in kN/m (=N/mm) (not
             functional at this timepoint!)
             False - free

        In its current iteration the supports may only be fixed or free for the
        given degree of freedom.
        """
        if self.supports is not None:
            ID = len(self.supports)+1
            support = Support(ID, node, t_x=t_x, t_y=t_y, t_z=t_z,
                              r_x=r_x, r_y=r_y, r_z=r_z)
            self.supports.append(support)
        else:
            ID = 1
            support = Support(ID, node, t_x=t_x, t_y=t_y, t_z=t_z,
                              r_x=r_x, r_y=r_y, r_z=r_z)
            self.supports = [support]

    def add_square_load(self, value, start_point, end_point,
                        direction=[0, 0, -1], beams=None, name='sq load',
                        auto=None):
        """
        Distributes a constant distributed load over a rectangular area over
        the beams forming the boundary of the area.

        value: Constant load value in N/mm
        start_point: Point class indicating the position of the first corner
        end_point: Point class indicatin the position of the second corner
        direction: Vector indicating the load direction
        beams: Edge beams of the area the load should be distributed to
        name: Name of the load
        auto: Indicates if the load should also be automatically distributed on
              all internal beams in the load area
        """

        sq_load = SquareLoad(value, start_point, end_point,
                             direction=direction, beams=beams, name=name)

        if auto is not None:
            sq_load.process_beams2(self)
        else:

            s_o = sq_load.process_beams(self)
            sq_load.distribute_load(s_o)

        if self.square_loads is not None:
            self.square_loads.append(sq_load)
        else:
            self.square_loads = [sq_load]

    def add_q_load(self, start_val, start_point, end_point, name=None,
                   end_val=None, direction=[0, 0, -1], beam=None):
        """
        Adds a distributed load on elements between two points.
        Load directions are added in the global coordinate system.
        start_val: starting value of the distributed load in N/mm
        end_val: end value of the distributed load in N/mm (optional, only to
                 be used if the load varies linearly)
        start_point: Point class indicating the start point of the load
        end_point: Point class indicating the end point of the load
        direction: Vector indicating the load direction in the global cs
        beam: Beam class the load should be limited to
        """

        if name is None:
            if self.q_loads is not None:
                name = 'q_load ' + str(len(self.q_loads)+1)
            else:
                name = 'q_load 1'

        q_load = QLoad(start_val, start_point, end_point, name=name,
                       end_val=end_val, direction=direction, beam=beam)
        q_load.process_nodes(self)

        if self.q_loads is not None:
            self.q_loads.append(q_load)
        else:
            self.q_loads = [q_load]

        return q_load

    def add_p_load(self, node, name=None, Fx=0, Fy=0, Fz=0, Tx=0,
                   Ty=0, Tz=0):
        """
        Adds a point load at a certain node.
        Load directions are added in the global coordinate system.

        node: Node class the load should be applied to
        name: Name of the load
        Fx-Tz: Forces and moments in/around the global cs directions
        """

        if self.p_loads is not None:
            if name is None:
                ID = len(self.p_loads)+1
                load = PointLoad(node, ID, Fx, Fy, Fz, Tx, Ty, Tz)
                load.apply()
            else:
                load = PointLoad(node, name, Fx, Fy, Fz, Tx, Ty, Tz)
                load.apply()
            self.p_loads.append(load)
        else:
            if name is None:
                ID = 1
                load = PointLoad(node, ID, Fx, Fy, Fz, Tx, Ty, Tz)
                load.apply()
            else:
                load = PointLoad(node, name, Fx, Fy, Fz, Tx, Ty, Tz)
                load.apply()
            self.p_loads = [load]
        return load

    def clear_weight(self):
        """
        Clear all the self-weight loads from the system.
        """

        sw = [l for l in self.q_loads if 'Self weight' in l.name]
        for q_load in sw:
            for p_load in q_load.loads:
                p_load.remove()
