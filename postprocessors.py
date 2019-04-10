import numpy as np


class PostProcessor():
    """
    The PostProcessor class is used after the structural system has been solved
    and performs the actions necessary to determine the results of the
    calculation.

    In its current iteration, only the BeamProcessor is available.
    """

    def __init__(self, system):
        """
        Initializes the PostProcessor, assigning the System class it operates
        in.
        """
        self.system = system


class BeamPostProcessor(PostProcessor):
    """
    The BeamPostProcessor is the postprocessor used when there are strictly
    Beam elements used in the system.
    """

    def __init__(self, system):
        """
        Initializes the BeamPostProcessor, assigning the System class it
        operates in.
        """
        PostProcessor.__init__(self, system)

    def assign_displacements(self):
        """
        Assigns the calculated displacements to the local nodes.
        """
        try:
            self.system.u
        except AttributeError:
            raise Exception("The system hasn't been solved yet, why is the" +
                            " postprocessor assigning displacements?")
            return
        for node in self.system.nodes:
            n = (node.ID - 1)*6
            node.displacement_results = np.array([self.system.u[n],
                                                  self.system.u[n+1],
                                                  self.system.u[n+2],
                                                  self.system.u[n+3],
                                                  self.system.u[n+4],
                                                  self.system.u[n+5]])

        for beam in self.system.beams:
            beam.displacements()

    def assign_forces(self, s=6):
        """
        Calculates and assigns the internal forces in the system to the local
        nodes.
        """
        for element in self.system.elements:
            element.internal_forces(s=s)
        for beam in self.system.beams:
            beam.moments(s=s)
            beam.normal_forces(s=s)
            beam.shear_forces(s=s)
