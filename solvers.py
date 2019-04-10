import numpy as np
import scipy


class Solver():
    """
    The Solver class is the parent class of (in the future) different solvers
    to solve the systems of equations assembled by the different Assembler
    classes.

    In its current iteration only the SimpleSolver is available.
    """

    def __init__(self, system):
        """
        Initiates the Solver class and assigns the parent System class using
        the system argument.
        """
        self.system = system


class SimpleSolver(Solver):
    """
    The SimpleSolver class.
    """

    def __init__(self, system):
        """
        Initiates the SimpleSolver class and assigns the parent System class
        using the system argument.
        """
        Solver.__init__(self, system)

    def solve(self, t_test=False):
        """
        This method solves the system of equations contained in the System
        class specified.

        It is possible test the time taken to solve the system of equations by
        setting the t_test argument to True.
        """

        print()
        print("Solving system of equations")
        print()

        if t_test is True:
            import time
            t_t = time.perf_counter()

        k_r = self.system.k_reduced
        f_r = self.system.f_reduced

        sol_ = scipy.linalg.solve(k_r, f_r, assume_a='pos')

        for i_ in reversed(range(33)):

            if i_ >= 16:

                tol = 10**-(i_ - 16)

            else:

                tol = 10**(10**(16 - i_))

            if np.allclose(np.dot(k_r, sol_), f_r, atol=tol):

                print()
                print("    Solution accuracy {:.1e}".format(tol))
                print()

                self.system.u_red = sol_

                sol = np.zeros((self.system.k_matrix.shape[0],))
                c = 0
                for i in self.system.assembler.remain_indices:
                    sol[i] = sol_[c]
                    c += 1

                self.system.u = sol
                f = np.dot(self.system.k_matrix, self.system.u)

                self.system.f = f
                if t_test is True:
                    t_t = time.perf_counter() - t_t
                    print()
                    print("        Solving with {:.1e}".format(tol) +
                          " accuracy took: {:.2f} s".format(t_t))
                    print()
                return
