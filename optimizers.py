import numpy as np
from edsFEM.profiles import I_Beams
from edsFEM.sections import IBeam


class Optimizer():
    """
    The Optimizers will in the future allow for automatic optimization of
    structural components.
    """

    def __init__(self, system):
        """
        Initializes the Optimizer class, assigning the System class it operates
        in.
        """
        self.sys = system


class BeamOptimizer(Optimizer):
    """
    The BeamOptimizer class optimizes beams used in the given system based
    on a very limited simplified Von Mises stress check.
    It currently only operates using the standard profile library in the
    edsFEM.profiles module and on systems using these standard profiles.
    """

    def __init__(self, system):
        """
        Initializes the BeamOptimizer class, assigning the System class it
        operates in.
        """
        Optimizer.__init__(self, system)

    def standard_sections(self, s_max=235):
        """
        This method checks if it is possible to replace the beam profiles in
        the system with a profile weighing less. If it is possible, the beams'
        profile is replaced with the lighter profile.
        """

        for beam in self.sys.beams:
            if beam.profile is not None:
                if type(beam.profile) is IBeam:
                    v_mis = 0

                    for el in beam.elements:
                        try:
                            tau = np.max(el.tau)
                        except AttributeError:
                            tau = 0
                        try:
                            sigma = np.max(el.sigma)
                        except AttributeError:
                            sigma = 0

                        if sigma + (tau / np.sqrt(3)) > v_mis:
                            v_mis = sigma + (tau / np.sqrt(3))
                            el_max = el

                    bigger_profs = [p for p in I_Beams if
                                    self.check_prof(p, el_max, s_max=s_max)
                                    ]

                    bigger_profs = sorted(bigger_profs, key=lambda x: x.g)

                    if bigger_profs:
                        beam.set_profile(bigger_profs[0])
                    else:
                        new_prof = sorted(I_Beams, key=lambda x: x.g)
                        new_prof = new_prof[-1]
                        beam.set_profile(new_prof)
                        print('better prof not found')

    def check_prof(self, profile, element, s_max=235):
        """
        This method checks if the simplified Von Mises stress check is satis-
        fied for the given element and profile.
        """
        p = profile

        S_z = (0.5 * (p.h-2*p.t_f) * p.t_w) * 0.25 * (p.h-2*p.t_f)
        S_z += p.b * p.t_f * (0.5 * p.h - 0.5 * p.t_f)
        S_y = ((p.b-0.5*p.t_w) * p.t_f) * 0.25 * (p.b)
        S_y += (p.h-2*p.t_f)*0.5*p.t_w*0.25*p.t_w

        tau_mz = (np.max(np.abs(element.Vz)) * S_z) / (p.I_y * p.t_w)
        tau_my = (np.max(np.abs(element.Vy)) * S_y) / (p.I_z * p.t_w)

        tau = tau_mz + tau_my

        sigma_normal = np.max(np.abs(element.N) / p.A)
        m_y = np.max(np.abs(element.My) / p.Wy_el)
        m_z = np.max(np.abs(element.Mz) / p.Wz_el)
        sigma_bending = np.sum([m_y, m_z])

        sigma = sigma_normal + sigma_bending

        v_mis = sigma + (tau / np.sqrt(3))

        if v_mis < s_max:
            return True
        else:
            return False
