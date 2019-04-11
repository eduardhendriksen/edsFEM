import numpy as np


class Section():
    """
    The section class stores the basic mechanical properties of a structural
    section.
    """

    def __init__(self, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el):
        """
        Initializes the Section class.

        E: Young's modulus (N/mm²)
        A: Area (mm²)
        I_y: Moment of Inertia around y-axis (mm^4)
        I_z: Moment of Inertia around z-axis (mm^4)
        G: Shear modulus (N/mm²)
        I_t: torsion property of the section (mm^4)
        Wy_el: Elastic section modulus about the y-axis (mm^3)
        Wz_el: Elastic section modulus about the z-axis (mm^3)
        """
        self.E = E
        self.A = A
        self.I_y = I_y
        self.I_z = I_z
        self.G = G
        self.I_t = I_t
        self.Wy_el = Wy_el
        self.Wz_el = Wz_el
        self.g_stiff = 0
        self.stiffeners = None
        self.type = None


class IBeam(Section):
    """
    The IBeam class stores the mechanical properties of the I beam type of
    structural section. It also contains methods pertaining to this type of
    steel member.
    """

    def __init__(self, E, G, designation, g, h, b, t_w, t_f, r, A, I_y, Wy_el,
                 Wy_pl, i_y, A_vz, I_z, Wz_el, Wz_pl, i_z, S_s, I_t, I_w,
                 stiffeners=None):
        """
        Initializes the IBeam class.

        E: Young's modulus (N/mm²)
        G: Shear modulus (N/mm²)
        designation: Name of the section
        g: Weight of the section (kg/mm)
        h: Section height (mm)
        b: Section width (mm)
        t_w: Web thickness (mm)
        t_f: Flange thickness (mm)
        r: Radius of root fillet, or toe radius (mm)
        A: Area (mm²)
        I_y: Moment of Inertia around y-axis (mm^4)
        Wy_el: Elastic section modulus about the y-axis (mm^3)
        Wy_pl: Plastic section modulus about the y-axis (mm^3)
        i_y: Radius of giration about y-axis (mm)
        A_vz: Shear area for shear in z-direction (mm²)
        I_z: Moment of Inertia around z-axis (mm^4)
        Wz_el: Elastic section modulus about the z-axis (mm^3)
        Wz_pl: Plastic section modulus about the z-axis (mm^3)
        i_z: Radius of giration about z-axis (mm)
        S_s: Statical moment of area about z-axis (mm^3)
        I_t: Torsion property of the section (mm^4)
        I_w: Warping constant (mm^6)
        """

        Section.__init__(self, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el)
        self.designation = designation
        self.g = g
        self.h = h
        self.b = b
        self.t_w = t_w
        self.t_f = t_f
        self.r = r
        self.Wy_pl = Wy_pl
        self.i_y = i_y
        self.A_vz = A_vz
        self.Wz_pl = Wz_pl
        self.i_z = i_z
        self.S_s = S_s
        self.I_w = I_w
        self.type = 'IBeam'

        self.stiffeners = stiffeners

        if self.stiffeners is not None:
            self.add_stiffeners()

    def __repr__(self):
        return ("%s\n" % (self.designation))

    def add_stiffeners(self, stiffeners=None):
        """
        Method to add a Stiffeners class to the section using the stiffeners
        argument.
        The relevant sectional properties are recalculated with the stiffeners
        applied.
        """
        if stiffeners is None:
            if self.stiffeners is not None:
                if type(stiffeners) is not list:
                    stiffeners = [stiffeners]
                for stiff in self.stiffeners:
                    if stiff.direction == 'hor':
                        th = stiff.thickness
                        b = self.b + 2 * th
                        h = self.h
                        self.I_y += (((1 / 12) * (b) * h**3) -
                                     ((1 / 12) * (self.b) * h**3))
                        self.Wy_el += ((((1 / 12) * (b) * h**3) -
                                        ((1 / 12) * (self.b) * h**3)) / (h/2))
                        self.I_z += (((1 / 12) * (h) * b**3) -
                                     ((1 / 12) * (h) * (self.b)**3))
                        self.Wz_el += ((((1 / 12) * (h) * b**3) -
                                       ((1 / 12) * (h) * (self.b)**3)) / (b/2))
                        self.g_stiff = 2 * (h) * (th) * 7800e-9
                        self.g += self.g_stiff

                    if stiff.direction == 'hor':
                        self.designation += (' s. ' +
                                             str(int(stiff.thickness)) +
                                             'mm long.')
                    else:
                        raise NotImplementedError("Only horizontal stiffeners "
                                                  + "implemented")
        else:
            if type(stiffeners) is not list:
                stiffeners = [stiffeners]
            self.stiffeners = stiffeners
            for stiff in stiffeners:
                if stiff.direction == 'hor':
                    th = stiff.thickness
                    b = self.b + 2 * th
                    h = self.h
                    self.I_y += (((1 / 12) * (b) * h**3) -
                                 ((1 / 12) * (self.b) * h**3))
                    self.Wy_el += ((((1 / 12) * (b) * h**3) -
                                    ((1 / 12) * (self.b) * h**3)) / (h/2))
                    self.I_z += (((1 / 12) * (h) * b**3) -
                                 ((1 / 12) * (h) * (self.b)**3))
                    self.Wz_el += ((((1 / 12) * (h) * b**3) -
                                   ((1 / 12) * (h) * (self.b)**3)) / (b/2))
                    self.g_stiff = 2 * (h) * (th) * 7800e-9
                    self.g += self.g_stiff
                else:
                    raise NotImplementedError("Only horizontal stiffeners "
                                              + "implemented")


class WeldedIBeam(IBeam):
    """
    The class containing the sectional properties and methods to determine them
    for built-up welded IBeam sections.
    """

    def __init__(self, E, G, h, b, t_w, t_f, a, stiffeners=None,
                 designation=None):
        """
        Initializes the WeldedIBeam class.

        E: Young's modulus (N/mm²)
        G: Shear modulus (N/mm²)
        h: Section height (mm)
        b: Section width (mm)
        t_w: Web thickness (mm)
        t_f: Flange thickness (mm)
        a: Throat thickness of the fillet welds (mm)
        stiffeners: Stiffeners class of stiffeners to be added to the section
        designation: Name of the section
        """
        r = a
        A = self.calc_A(h, b, t_w, t_f, a)
        g = self.calc_g(A)
        I_y, I_z = self.calc_I(h, b, t_w, t_f, a)
        Wy_el, Wy_pl, Wz_el, Wz_pl = self.calc_W(h, b, t_w, t_f, a)
        i_y, i_z = self.calc_i(h, b, t_w, t_f, a)
        A_vz = self.calc_Avz(h, t_w)
        S_s = self.calc_S_s(h, b, t_w, t_f, a)
        I_t = self.calc_I_t(h, b, t_w, t_f, a)
        I_w = self.calc_I_w(h, b, t_w, t_f, a)
        if designation is None:
            designation = ('Welded: ' + str(round(h, 1)) + 'x' +
                           str(round(b, 1)) +
                           ' t_f ' + str(round(t_f, 1)) + ' t_w ' +
                           str(round(t_w, 1))
                           )
        IBeam.__init__(self, E, G, designation, g, h, b, t_w, t_f, r, A, I_y,
                       Wy_el, Wy_pl, i_y, A_vz, I_z, Wz_el, Wz_pl, i_z,
                       S_s, I_t, I_w, stiffeners=stiffeners)

    def calc_A(self, h, b, t_w, t_f, a):
        """
        Calculates the area of the welded section in mm² based on the arguments

        h: Section height (mm)
        b: Section width (mm)
        t_w: Web thickness (mm)
        t_f: Flange thickness (mm)
        a: Throat thickness of the fillet welds (mm)

        returns:

        A: Area of the welded section (mm²)
        """
        A = 2 * b * t_f
        A += (h - 2 * t_f) * t_w
        A += 4 * (a**2 / 2)

        return A

    def calc_g(self, A):
        """
        Calculates the weight of the welded section in kg/mm based on the Area
        specified in the A argument (mm²)

        A: Cross-sectional area of the section (mm²)

        returns:

        g: Weight of the section (kg/mm)
        """

        rho = 7800e-9

        return A * rho

    def calc_I(self, h, b, t_w, t_f, a):
        """
        Calculates I_y and I_z based on the values specified in the arguments.

        h: Section height (mm)
        b: Section width (mm)
        t_w: Web thickness (mm)
        t_f: Flange thickness (mm)
        a: Throat thickness of the fillet welds (mm)

        returns:

        I_y: Moment of Inertia around y-axis (mm^4)
        I_z: Moment of Inertia around y-axis (mm^4)
        """

        I_y = ((1/12) * b * h**3) - ((1/12) * (b-t_w) * (h-2*t_f)**3)
        I_z = ((1/12) * h * b**3) - ((1/12) * (h-2*t_f) * (b-t_w)**3)

        return I_y, I_z

    def calc_W(self, h, b, t_w, t_f, a):
        """
        Calculates Wy_el, Wy_pl, Wz_el and Wz_pl based on the values specified
        in the arguments.

        h: Section height (mm)
        b: Section width (mm)
        t_w: Web thickness (mm)
        t_f: Flange thickness (mm)
        a: Throat thickness of the fillet welds (mm)

        returns:

        Wy_el: Elastic section modulus about the y-axis (mm^3)
        Wy_pl: Plastic section modulus about the y-axis (mm^3)
        Wz_el: Elastic section modulus about the z-axis (mm^3)
        Wz_pl: Plastic section modulus about the z-axis (mm^3)
        """

        I_y = ((1/12) * b * h**3) - ((1/12) * (b-t_w) * (h-2*t_f)**3)
        Wy_el = I_y / (h / 2)
        Wy_pl = 2 * ((b*t_f * (h/2-t_f/2)) + ((h/2-t_f)*t_w * ((h/2-t_f)/2)))
        I_z = ((1/12) * h * b**3) - ((1/12) * (h-2*t_f) * (b-t_w)**3)
        Wz_el = I_z / (b / 2)
        Wz_pl = 2 * (((b/2)*t_f * (b/4)) + ((h-2*t_f)*(t_w/2) * (t_w/4)))

        return Wy_el, Wy_pl, Wz_el, Wz_pl

    def calc_i(self, h, b, t_w, t_f, a):
        """
        Calculates i_y and i_z based on the values specified in the arguments.

        h: Section height (mm)
        b: Section width (mm)
        t_w: Web thickness (mm)
        t_f: Flange thickness (mm)
        a: Throat thickness of the fillet welds (mm)

        returns:

        i_y: Radius of giration about y-axis (mm)
        i_z: Radius of giration about z-axis (mm)
        """

        A = 2 * b * t_f + (h - 2 * t_f) * t_w + 4 * (a**2 / 2)
        I_y = ((1/12) * b * h**3) - ((1/12) * (b-t_w) * (h-2*t_f)**3)
        i_y = I_y / A
        I_z = ((1/12) * h * b**3) - ((1/12) * (h-2*t_f) * (b-t_w)**3)
        i_z = I_z / A

        return i_y, i_z

    def calc_Avz(self, h, t_w):
        """
        Calculates A_vz based on the values specified in the arguments.
        Based on EN 1993-1-5 DE, hochbau.

        h: Section height (mm)
        t_w: Web thickness (mm)

        returns:

        A_vz: Shear area for shear in z-direction (mm²)
        """

        eta = 1.2
        Avz = eta * h * t_w

        return Avz

    def calc_S_s(self, h, b, t_w, t_f, a):
        """
        Calculates S_s based on the values specified in the arguments.

        Does not function as yet!

        h: Section height (mm)
        b: Section width (mm)
        t_w: Web thickness (mm)
        t_f: Flange thickness (mm)
        a: Throat thickness of the fillet welds (mm)

        returns:

        S_s: Statical moment of area about z-axis (mm^3)
        """
        return 0

    def calc_I_t(self, h, b, t_w, t_f, a):
        """
        Calculates the torsional constant I_t

        h: Section height (mm)
        b: Section width (mm)
        t_w: Web thickness (mm)
        t_f: Flange thickness (mm)
        a: Throat thickness of the fillet welds (mm)

        returns:

        I_t: Torsional constant (mm^4)
        """
        I_t = 1.15 * (1/3) * ((2 * b * t_f**3) + ((h - 2*t_f)*t_w**3))
        return I_t

    def calc_I_w(self, h, b, t_w, t_f, a):
        """
        Calculates I_w based on the values specified in the arguments.

        Does not function as yet!

        h: Section height (mm)
        b: Section width (mm)
        t_w: Web thickness (mm)
        t_f: Flange thickness (mm)
        a: Throat thickness of the fillet welds (mm)

        returns:

        I_w: Warping constant (mm^6)
        """

        I_w = (1/24) * t_f * b**3 * (h-t_f)**2
        return I_w


class UBeam(Section):
    """
    The IBeam class stores the mechanical properties of the U beam type of
    structural section.
    """

    def __init__(self, E, G, designation, g, h, b, s, t, r, A, I_y, Wy_el,
                 Wy_pl, i_y, A_vz, I_z, Wz_el, Wz_pl, i_z, S_s, I_t, I_w):
        """
        Initializes the UBeam class.

        E: Young's modulus (N/mm²)
        G: Shear modulus (N/mm²)
        designation: Name of the section
        g: Weight of the section (kg/mm)
        h: Section height (mm)
        b: Section width (mm)
        s: Web thickness (mm)
        t: Flange thickness (mm)
        r: Radius of root fillet, or toe radius (mm)
        A: Area (mm²)
        I_y: Moment of Inertia around y-axis (mm^4)
        Wy_el: Elastic section modulus about the y-axis (mm^3)
        Wy_pl: Plastic section modulus about the y-axis (mm^3)
        i_y: Radius of giration about y-axis (mm)
        A_vz: Shear area for shear in z-direction (mm²)
        I_z: Moment of Inertia around z-axis (mm^4)
        Wz_el: Elastic section modulus about the z-axis (mm^3)
        Wz_pl: Plastic section modulus about the z-axis (mm^3)
        i_z: Radius of giration about z-axis (mm)
        S_s: Statical moment of area about z-axis (mm^3)
        I_t: Torsion property of the section (mm^4)
        I_w: Warping constant (mm^6)
        """

        Section.__init__(self, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el)
        self.designation = designation
        self.g = g
        self.h = h
        self.b = b
        self.s = s
        self.t = t
        self.r = r
        self.Wy_pl = Wy_pl
        self.i_y = i_y
        self.A_vz = A_vz
        self.Wz_pl = Wz_pl
        self.i_z = i_z
        self.S_s = S_s
        self.I_w = I_w
        self.type = 'UBeam'


class LBeam(Section):
    """
    The LBeam class stores the mechanical properties of the L beam type of
    structural section.
    """

    def __init__(self, E, G, designation, g, h, b, t, r_2, r_1, A, I_y,
                 Wy_el, i_y, I_z, Wz_el, i_z, I_t):
        """
        Initializes the LBeam class.

        E: Young's modulus (N/mm²)
        G: Shear modulus (N/mm²)
        designation: Name of the section
        g: Weight of the section (kg/mm)
        h: Section height (mm)
        b: Section width (mm)
        t: Flange thickness (mm)
        r_2: Toe radius (mm)
        r_1: Radius of root fillet (mm)
        A: Area (mm²)
        I_y: Moment of Inertia around y-axis (mm^4)
        Wy_el: Elastic section modulus about the y-axis (mm^3)
        i_y: Radius of giration about y-axis (mm)
        I_z: Moment of Inertia around z-axis (mm^4)
        Wz_el: Elastic section modulus about the z-axis (mm^3)
        i_z: Radius of giration about z-axis (mm)
        I_t: Torsion property of the section (mm^4)
        """

        Section.__init__(self, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el)
        self.designation = designation
        self.g = g
        self.h = h
        self.b = b
        self.t = t
        self.r_2 = r_2
        self.r_1 = r_1
        self.i_y = i_y
        self.i_z = i_z
        self.type = 'LBeam'


class Stiffeners():
    """
    The Stiffeners class contains mechanical properties of stiffeners to be
    applied to the different sections that allow them.
    """

    def __init__(self, direction, thickness):
        """
        Initializes the Stiffeners class.

        direction: 'hor' for horizontally applied stiffeners (longitudinal)
                   Any other direction is not currently available.
        thickness: Thickness of the stiffener plates (mm)
        """

        self.direction = direction
        self.thickness = thickness


class CHS(Section):
    """
    The CHS class stores the mechanical properties and methods of the Circular
    Hollow Section type of structural section.
    """

    def __init__(self, E, G, D, t, rho=7800e-9):
        """
        Initializes the CHS class and calculates the relevant mechanical
        properties of the section.

        E: Young's modulus (N/mm²)
        G: Shear modulus (N/mm²)
        D: Section diameter (mm)
        t: Section thickness (mm)
        rho: Section density (kg/mm^3)
        """
        self.designation = "Ø {:.1f} mm, t {:.1f} mm".format(D, t)
        self.D = D
        self.t = t
        A = self.calc_A(D, t)
        self.g = rho * A
        I_y, I_z = self.calc_I(D, t)
        I_t = self.calc_I_t(D, t)
        Wy_el, Wz_el = self.calc_W(D, t)
        self.type = 'CHS'
        Section.__init__(self, E, A, I_y, I_z, G, I_t, Wy_el, Wz_el)

    def calc_A(self, D, t):
        """
        Calculates the area of the CHS

        D: Section diameter (mm)
        t: Section thickness (mm)

        returns:
        A: Area of the CHS (mm²)
        """
        return (1/4)*np.pi*(D)**2 - (1/4)*np.pi*(D-2*t)**2

    def calc_I(self, D, t):
        """
        Calculates the moments of inertia of the CHS

        D: Section diameter (mm)
        t: Section thickness (mm)

        returns:
        I_y: Moment of Inertia around y-axis (mm^4)
        I_z: Moment of Inertia around y-axis (mm^4)
        """
        I_y_z = (1/64)*np.pi*D**4 - (1/64)*np.pi*(D-2*t)**4
        return I_y_z, I_y_z

    def calc_W(self, D, t):
        """
        Calculates the elastic section moduli of the CHS

        D: Section diameter (mm)
        t: Section thickness (mm)

        returns:
        Wy_el: Elastic section modulus around y-axis (mm^4)
        Wz_el: Elastic section modulus around y-axis (mm^4)
        """
        I_y_z = (1/64)*np.pi*D**4 - (1/64)*np.pi*(D-2*t)**4
        W_el = (2 * I_y_z) / (D-2*t)
        return W_el, W_el

    def calc_I_t(self, D, t):
        """
        Calculates the torsional constant of the CHS

        D: Section diameter (mm)
        t: Section thickness (mm)

        returns:
        I_t: Torsional constant (mm^4)
        """
        return 2 * ((1/64)*np.pi*D**4 - (1/64)*np.pi*(D-2*t)**4)
