import numpy as np


class StructureNode():
    """
    Special class of node, placed at the ends of beams.
    These nodes are rendered by the renderer.
    """

    def __init__(self, ID, node=None, point=None):
        """
        Initializes the StructureNode class.

        ID: integer representing the number of the StructureNode
        node: Node class the StructureNode is placed on
        point: Point class indicating the location of the StructureNode
        """

        self.ID = ID
        if node is None:
            if point is None:
                raise ValueError('Initiate a structure node at a location')
            else:
                self.point = point
        else:
            self.node = node
            self.point = self.node.point
            self.node.structurenode = self

    def __repr__(self):
        return("\nNode\n"
               "\nID = {}\n"
               "node ID = {}\n"
               "Point = {}\n"
               "Fx = {:.2f}\n"
               "Fy = {:.2f}\n"
               "Fz = {:.2f}\n"
               "Tx = {:.2f}\n"
               "Ty = {:.2f}\n"
               "Tz = {:.2f}\n"
               "ux = {:.2f}\n"
               "uy = {:.2f}\n"
               "uz = {:.2f}\n"
               "phi_x = {:.2f}\n"
               "phi_y = {:.2f}\n"
               "phi_z = {:.2f}\n".format(self.ID, self.node.ID, self.point,
                                         *self.node.force_results,
                                         *self.node.displacement_results))

    def __str__(self):
        return("\nNode\n"
               "\nID = {}\n"
               "node ID = {}\n"
               "Point = {}\n"
               "Fx = {:.2f}\n"
               "Fy = {:.2f}\n"
               "Fz = {:.2f}\n"
               "Tx = {:.2f}\n"
               "Ty = {:.2f}\n"
               "Tz = {:.2f}\n"
               "ux = {:.2f}\n"
               "uy = {:.2f}\n"
               "uz = {:.2f}\n"
               "phi_x = {:.2f}\n"
               "phi_y = {:.2f}\n"
               "phi_z = {:.2f}\n".format(self.ID, self.node.ID, self.point,
                                         *self.node.force_results,
                                         *self.node.displacement_results))


class Node:
    """
    The Node class represents the interface points between the elements out of
    which the structural system in built up.
    At Nodes, supports and concentrated loads can be applied and structural
    elements can be connected.
    """

    def __init__(self, ID, point, Fx=0, Fy=0, Fz=0, Tx=0, Ty=0, Tz=0,
                 ux=0, uy=0, uz=0, phi_x=0, phi_y=0, phi_z=0,
                 support=None, elements=None):
        """
        ID: ID of the node, integer
        point: Point object
        Fx: Value of Fx
        Fy: Value of Fy
        Fz: Value of Fz
        Tx: Value of Tx (Moment around x-axis)
        Ty: Value of Ty (Moment around z-axis)
        Tz: Value of Tz (Moment around y-axis)
        ux: Value of ux
        uy: Value of uy
        uz: Value of uz
        phi_x: Value of phi_x (rotation around x-axis)
        phi_y: Value of phi_y (rotation around z-axis)
        phi_z: Value of phi_z (rotation around y-axis)
        support: Support object
        elements: List of attached elements
        """
        self.ID = ID
        self.point = point

        self.Fx = Fx
        self.Fy = Fy
        self.Fz = Fz
        self.Tx = Tx
        self.Ty = Ty
        self.Tz = Tz

        self.ux = ux
        self.uy = uy
        self.uz = uz
        self.phi_x = phi_x
        self.phi_y = phi_y
        self.phi_z = phi_z

        self.support = support

        self.elements = elements

        self.force_results = np.zeros((6,))
        self.displacement_results = np.zeros((6,))

    def __repr__(self):
        return("\nNode\n"
               "\nID = {}\n"
               "Point = {}\n"
               "Fx = {:.2f}\n"
               "Fy = {:.2f}\n"
               "Fz = {:.2f}\n"
               "Tx = {:.2f}\n"
               "Ty = {:.2f}\n"
               "Tz = {:.2f}\n"
               "ux = {:.2f}\n"
               "uy = {:.2f}\n"
               "uz = {:.2f}\n"
               "phi_x = {:.2f}\n"
               "phi_y = {:.2f}\n"
               "phi_z = {:.2f}\n".format(self.ID, self.point,
                                         *self.force_results,
                                         *self.displacement_results))

    def __str__(self):
        return("\nNode\n"
               "\nID = {}\n"
               "Point = {}\n"
               "Fx = {:.2f}\n"
               "Fy = {:.2f}\n"
               "Fz = {:.2f}\n"
               "Tx = {:.2f}\n"
               "Ty = {:.2f}\n"
               "Tz = {:.2f}\n"
               "ux = {:.2f}\n"
               "uy = {:.2f}\n"
               "uz = {:.2f}\n"
               "phi_x = {:.2f}\n"
               "phi_y = {:.2f}\n"
               "phi_z = {:.2f}\n".format(self.ID, self.point,
                                         *self.force_results,
                                         *self.displacement_results))
