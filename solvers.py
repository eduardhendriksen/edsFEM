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

        k = self.system.k_matrix
        f = self.system.f_vector

        sol = scipy.linalg.solve(k, f, assume_a='pos')

        for i_ in reversed(range(33)):

            if i_ >= 16:

                tol = 10**-(i_ - 16)

            else:

                tol = 10**(10**(16 - i_))

            if np.allclose(np.dot(k, sol), f, atol=tol):

                print()
                print("    Solution accuracy {:.1e}".format(tol))
                print()

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
