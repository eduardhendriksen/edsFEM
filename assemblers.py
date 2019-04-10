import numpy as np


class Assembler():
    """
    The Assembler classes formulate the stiffness matrix and force vector
    for which the system of equations is solved.

    At this point in time, only the BeamAssembler class is available.
    """

    def __init__(self, system):
        """
        Initializes the Assembler class and assigns the System class calling
        the Assembler.
        """
        self.system = system
        self.removed_indices = []
        self.remain_indices = []


class BeamAssembler(Assembler):
    """
    The BeamAssembler class formulates the stiffness matrix and force vector
    for which the system of equations is solved for problems in which Beam
    elements are utilized.
    """

    def __init__(self, system):
        """
        Initializes the BeamAssembler class and assigns the System class
        calling the BeamAssembler.
        """
        Assembler.__init__(self, system)

    def assemble(self, selfweight=None, t_test=False):
        """
        Assembles a symmetric stiffness matrix and removes restrained degrees
        of freedom from the equilibrium equations.

        Also assigns self-weight using the self.weight method.
        """
        print()
        print("Assembling stiffness matrix")
        print()
        if t_test is True:

            import time

        if selfweight is not None:

            if t_test is False:

                self.weight(selfweight)

            else:

                t_w = time.perf_counter()
                self.weight(selfweight)

                t_w = time.perf_counter() - t_w

                print()
                print("    Assigning member weight took: " +
                      "{:.2f}".format(t_w) + " s")
                print()

        if t_test is True:

            t_r = time.perf_counter()

        for element in self.system.elements:

            element.releases()

        if t_test is True:

            t_r = time.perf_counter() - t_r

            print()
            print("    Assigning member releases took: " +
                  "{:.2f}".format(t_r) + " s")
            print()

        if t_test is True:

            t_b = time.perf_counter()

        if self.system.beams is not None:

            for beam in self.system.beams:

                beam.sort_elements()

        if t_test is True:

            t_b = time.perf_counter() - t_b

            print()
            print("    Sorting beam elements took: " +
                  "{:.2f}".format(t_b) + " s")
            print()

        if t_test is True:

            t_a = time.perf_counter()

        dofs = self.system.dims*2
        shape = len(self.system.nodes) * dofs

        k = np.zeros((shape, shape))
        f = np.zeros((shape,))
        u = np.zeros((shape,))

        for node, i in zip(self.system.nodes,
                           list(range(len(self.system.nodes)))):

            for element in node.elements:

                k_ = element.glob_stiffness_matrix

                if node is element.node_1:

                    j_1 = (node.ID-1)*dofs
                    j_2 = (element.node_2.ID-1)*dofs

                    for m in range(dofs):

                        for n in range(dofs):

                            k[j_1+m, j_1+n] += k_[m, n]

                    for m in range(dofs):

                        for n in range(dofs):

                            k[j_1+m, j_2+n] += k_[m, dofs+n]

                    for m in range(dofs):

                        for n in range(dofs):

                            k[j_2+m, j_1+n] += k_[dofs+m, n]

                    for m in range(dofs):

                        for n in range(dofs):

                            k[j_2+m, j_2+n] += k_[dofs+m, dofs+n]

            j = (node.ID-1) * dofs

            if node.support is not None:

                if node.support.t_x is 0:

                    u[j] = 0

                    if node.ux != 0:

                        u[j] += node.ux

                    f[j] = node.Fx

                elif node.support.t_x is False:

                    if node.ux is 0:

                        u[j] = np.NaN

                    else:

                        u[j] = node.ux

                    f[j] = node.Fx

                else:

                    if node.ux is 0:

                        u[j] = np.NaN

                    else:

                        u[j] = node.ux

                    k[j, j] += node.support.t_x
                    f[j] = node.Fx

                if node.support.t_y is 0:

                    u[j+1] = 0

                    if node.uy != 0:

                        u[j+1] += node.uy

                    f[j+1] = node.Fy

                elif node.support.t_y is False:

                    if node.uy is 0:

                        u[j+1] = np.NaN

                    else:

                        u[j+1] = node.uy

                    f[j+1] = node.Fy

                else:

                    if node.uy is 0:

                        u[j+1] = np.NaN

                    else:

                        u[j+1] = node.uy

                    k[j+1, j+1] += node.support.t_y
                    f[j+1] = node.Fy

                if node.support.t_z is 0:

                    u[j+2] = 0

                    if node.uz != 0:

                        u[j+2] += node.uz

                    f[j+2] = node.Fz

                elif node.support.t_z is False:

                    if node.uz is 0:

                        u[j+2] = np.NaN

                    else:

                        u[j+2] = node.uz

                    f[j+2] = node.Fz

                else:

                    if node.uz is 0:

                        u[j+2] = np.NaN

                    else:

                        u[j+2] = node.uz

                    k[j+2, j+2] += node.support.t_z
                    f[j+2] = node.Fz

                if node.support.r_x is 0:

                    u[j+3] = 0

                    if node.phi_x != 0:

                        u[j+3] += node.phi_x

                    f[j+3] = node.Tx

                elif node.support.r_x is False:

                    if node.phi_x is 0:

                        u[j+3] = np.NaN

                    else:

                        u[j+3] = node.phi_x

                    f[j+3] = node.Tx

                else:

                    if node.phi_x is 0:

                        u[j+3] = np.NaN

                    else:

                        u[j+3] = node.phi_x

                    k[j+3, j+3] += node.support.r_x
                    f[j+3] = node.Tx

                if node.support.r_y is 0:

                    u[j+4] = 0

                    if node.phi_y != 0:

                        u[j+4] += node.phi_y

                    f[j+4] = node.Ty

                elif node.support.r_y is False:

                    if node.phi_y is 0:

                        u[j+4] = np.NaN

                    else:

                        u[j+4] = node.phi_y

                    f[j+4] = node.Ty

                else:

                    if node.phi_y is 0:

                        u[j+4] = np.NaN

                    else:

                        u[j+4] = node.phi_y

                    k[j+4, j+4] += node.support.r_y
                    f[j+4] = node.Ty

                if node.support.r_z is 0:

                    u[j+5] = 0

                    if node.phi_z != 0:

                        u[j+5] += node.phi_z

                    f[j+5] = node.Tz

                elif node.support.r_z is False:

                    if node.phi_z is 0:

                        u[j+5] = np.NaN

                    else:

                        u[j+5] = node.phi_z

                    f[j+5] = node.Tz

                else:

                    if node.phi_z is 0:

                        u[j+5] = np.NaN

                    else:

                        u[j+5] = node.phi_z

                    k[j+5, j+5] += node.support.r_z
                    f[j+5] = node.Tz

            else:

                if node.ux is 0:

                    u[j] = np.NaN

                else:

                    u[j] = node.ux

                if node.uy is 0:

                    u[j+1] = np.NaN

                else:

                    u[j+1] = node.uy

                if node.uz is 0:

                    u[j+2] = np.NaN

                else:

                    u[j+2] = node.uz

                if node.phi_x is 0:

                    u[j+3] = np.NaN

                else:

                    u[j+3] = node.phi_x

                if node.phi_y is 0:

                    u[j+4] = np.NaN

                else:

                    u[j+4] = node.phi_y

                if node.phi_z is 0:

                    u[j+5] = np.NaN

                else:

                    u[j+5] = node.phi_z

                f[j] = node.Fx
                f[j+1] = node.Fy
                f[j+2] = node.Fz
                f[j+3] = node.Tx
                f[j+4] = node.Ty
                f[j+5] = node.Tz

        if t_test is True:
            t_t = time.perf_counter()

        for i_ in reversed(range(33)):

            if i_ >= 16:
                tol = 10**-(i_ - 16)
            else:
                tol = 10**(10**(16 - i_))

            if np.allclose((k.transpose()), k, atol=tol):

                print()
                print("    Stiffness matrix is symmetric with accuracy: " +
                      "{:.1e}".format(tol))
                print()

                self.system.k_matrix = k
                self.system.f_vector = f
                self.system.u_vector = u

                # TODO: insert non-zero boundary conditions

                if t_test is True:

                    t_rem = time.perf_counter()

                z = np.where(self.system.u_vector == 0)
                i = np.arange(k.shape[0])
                i = np.delete(i, z, axis=0)

                k_reduced = np.mat(k[np.ix_(i, i)])
                u_reduced = np.delete(u, z, axis=0)
                f_reduced = np.delete(f, z, axis=0)

                self.system.k_reduced = k_reduced
                self.system.u_reduced = u_reduced
                self.system.f_reduced = f_reduced
                self.removed_indices = z
                self.remain_indices = i

                if t_test is True:

                    t_rem = time.perf_counter() - t_rem

                    print()
                    print("        Reducing the equation set took: " +
                          "{:.2f}".format(t_rem) + " s")
                    print()

                break

        if t_test is True:

            print()
            print("    Finding stiffness matrix tolerance took: " +
                  "{:.2f}".format(time.perf_counter() - t_t) + " s")
            print()

        if t_test is True:

            t_a = time.perf_counter() - t_a

            print()
            print("Stiffness matrix assembly took: " +
                  "{:.2f}".format(t_a) + " s")
            print()

    def weight(self, factor=1):
        """
        Assigns self-weight to the system.
        The factor given indicates the multiplication factor used for the
        application of self-weight.
        """

        for beam in self.system.beams:

            q_g = beam.g * 9.81 * factor

            self.system.add_q_load(q_g, beam.point_1, beam.point_2,
                                   name=('Self weight beam ' +
                                         str(beam.designation)), beam=beam)
