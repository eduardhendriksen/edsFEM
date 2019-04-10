import numpy as np
from edsFEM.points import Point
from itertools import combinations


class Rectangle():
    """
    The Rectangle class represents a rectangle spanned by four Points.
    It also contains various methods for geometrical checks and calculations.
    """

    def __init__(self, p_1, p_3, p_2=None, p_4=None, parent=None,
                 children=None):
        """
        Initializes the Rectangle class/

        p_1: Point class indicating the first point of the rectangle
        p_3: Point class indicating the third point of the rectangle
        p_2: Point class indicating the second point of the rectangle
        p_4: Point class indicating the fourth point of the rectangle
        parent: Rectangle class this Rectangle class is contained in
        children: Rectangles classes contained by this Rectangle
        """
        if p_1 != p_3:
            self.A = p_1
            self.C = p_3
            if p_2 is None or p_4 is None:
                self.calc_BD()
            else:
                self.B = p_2
                self.D = p_4

            self.calc_sides()
            self.rect = self.check()

            if self.rect:

                self.parent = parent
                self.children = children
                self.c_children = None
                self.calc_area()
                self.smallest = True

        else:

            raise ValueError("Can't create a rectangle with coinciding points")

    def __repr__(self):
        return ("A = %s\n"
                "B = %s\n"
                "C = %s\n"
                "D = %s\n" % (self.A, self.B, self.C, self.D))

    def __str__(self):
        return ("A = %s\n"
                "B = %s\n"
                "C = %s\n"
                "D = %s\n" % (self.A, self.B, self.C, self.D))

    def __eq__(self, other):
        if (self.A == other.A and (self.B == other.B and self.C == other.C
                                   and self.D == other.D)):
            return True
        else:
            return False

    def calc_BD(self):
        """
        Calculates points B and D from the position of points A and C
        """

        ps = self.C.v - self.A.v

        sides = ps[np.abs(ps) > 0]

        ps_1 = np.where(ps == sides[0], sides[0], 0)
        ps_2 = np.where(ps == sides[1], sides[1], 0)

        self.B = Point(*(self.A.v + ps_2))
        self.D = Point(*(self.A.v + ps_1))

    def calc_sides(self):
        """
        Calculates the vectors and length of the sides of the rectangle.
        """

        self.AB = self.B - self.A
        self.AB_l = self.AB.modulus()
        self.BC = self.C - self.A
        self.BC_l = self.BC.modulus()
        self.CD = self.D - self.C
        self.CD_l = self.CD.modulus()
        self.DA = self.A - self.D
        self.DA_l = self.DA.modulus()
        self.AC = self.C - self.A
        self.AC_l = self.AC.modulus()
        self.BD = self.D - self.B
        self.BD_l = self.BD.modulus()

        self.d_list = [self.AB_l, self.BC_l, self.CD_l, self.DA_l, self.AC_l,
                       self.BD_l]
        self.p_list = [self.A, self.B, self.C, self.D]

    def check(self, *points):
        """
        Checks if this is a rectangle, or if the given points form a rectangle.
        """

        if not points:

            z = set(self.d_list)

            if len(z) == 3:
                return True
            elif np.isclose(np.min(self.d_list) * np.sqrt(2),
                            np.max(self.d_list)):
                return True
            else:
                return False

        else:

            if any(self.A.is_rect(*comb) is True for comb in
                   combinations(*points, r=4)):
                return True
            else:
                return False

    def calc_area(self):
        """
        Calculates the area of this rectangle
        """
        self.Area = self.AB_l * self.BC_l

    def order(self, points):
        """
        Orders the given points so a rectangle is formed.
        """
        points = list(points)
        ps = [points[0]]
        points.remove(ps[0])

        while len(ps) < 4:
            for p in points:
                if np.count_nonzero(p.v - ps[-1].v) == 1:
                    ps.append(p)
                    points.remove(p)

        return ps

    def find_children(self, p_list):
        """
        Finds children rectangles in the point list.
        """

        self.Area_filled = self.Area

        self.children = []
        self.c_children = []

        p_list = [p for p in p_list if
                  p.in_rect(self.A, self.B, self.C, self.D)]

        for comb in combinations(p_list, r=4):

            if self.check(comb):

                if not all(p in self.p_list for p in comb):

                    comb = self.order(comb)

                    if (not any(p.inside_rect(*k.p_list) for
                                k in self.children + self.c_children
                                for p in comb) or
                        not any(p.inside_rect(*comb) for k in
                                self.children + self.c_children
                                for p in k.p_list)):

                            if not any(all(p.in_rect(*k.p_list)
                                       for p in comb)
                                       for k in self.children +
                                       self.c_children):

                                if not any(all(p == k for k in l.p_list for
                                               p in comb)
                                           for l in self.children +
                                           self.c_children):

                                    self.children.append(
                                                         Rectangle(comb[0],
                                                                   comb[2],
                                                                   comb[1],
                                                                   comb[3],
                                                                   parent=self)
                                                         )

                                    self.smallest = False

                                    if self.parent is not None:

                                        (self.parent.c_children
                                         ).append(self.children[-1])

                                    self.children[-1].find_children(p_list)

                            else:

                                n_rect = Rectangle(comb[0], comb[2],
                                                   comb[1], comb[3])

                                for rect in self.children:
                                    if all(p.in_rect(*k.p_list) for
                                           k in self.children +
                                           self.c_children for p in comb):
                                        if n_rect.Area < rect.Area:
                                            self.children.remove(rect)
                                            n_rect.parent = self

                                            self.children.append(n_rect)

                                            self.smallest = False

                                            if self.parent is not None:

                                                (self.parent.c_children
                                                 ).append(self.children[-1])

                                            (self.children[-1]
                                             ).find_children(p_list)

                    try:
                        self.Area_filled -= self.children[-1].Area
                    except AttributeError:
                        print(self.children[-1])

        if self.c_children:
            if self.parent:
                for child in self.c_children:
                    if child not in self.parent.c_children:
                        self.parent.c_children.append(child)
