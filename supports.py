import numpy as np


class Support():
    """
    The Support class provides the software equivalent of the structural
    mechanical supports.
    In its current iteration it only possible to model supports fully fixed
    or fully free in the different degrees of freedom.
    """

    def __init__(self, ID, node, t_x=None, t_y=None, t_z=None,
                 r_x=None, r_y=None, r_z=None, loc_cs=None):
        """
        node: the node the support is attached to
        t_x: 0 - fixed
             any number != 0 - support stiffness in kN/m (=N/mm) (not
             functional at this timepoint!)
             None - free

        In its current iteration the supports may only be fixed or free for the
        given degree of freedom.
        """
        self.ID = ID
        self.node = node
        self.node.support = self
        self.t_x = t_x
        self.t_y = t_y
        self.t_z = t_z
        self.r_x = r_x
        self.r_y = r_y
        self.r_z = r_z
        self.loc_cs = None
        self.v = np.array([t_x, t_y, t_z, r_x, r_y, r_z])

    def reactions(self):
        """
        Calculates the support reactions - Currently not performed in this
        class, but in the PostProcessor.
        """

        pass
