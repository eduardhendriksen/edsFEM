import numpy as np
from itertools import combinations


class Point:
    """
    The Point class contains the x, y and z coordinates of the specified point.
    It also contains various methods for geometrical checks and calculations.
    """

    def __init__(self, x, y, z):
        """
        Initializes the Point class, assigning its x, y and z positions.
        """
        self.x = x
        self.y = y
        self.z = z
        self.v = np.array([x, y, z])

    def __repr__(self):
        return ('Point(' + str(self.x) + ', ' + str(self.y) + ', ' +
                str(self.z) + ')')

    def __str__(self):
        return ('Point(' + str(self.x) + ', ' + str(self.y) + ', ' +
                str(self.z) + ')')

    def print_values(self):
        print(self.x, self.y, self.z)

    def modulus(self):
        return np.sqrt(np.sum(np.square(self.v)))

    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return Point(x, y, z)

    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z
        return Point(x, y, z)

    def __truediv__(self, other):
        if type(other) is Point:
            x = self.x / other.x
            y = self.y / other.y
            z = self.z / other.z
            return Point(x, y, z)
        else:
            x = self.x / other
            y = self.y / other
            z = self.z / other
            return Point(x, y, z)

    def __eq__(self, other):
        if self.x == other.x and self.y == other.y and self.z == other.z:
            return True
        else:
            return False

    def on_line(self, p_s, p_e, length=False):
        """
        Checks if this point is on the line between given points p_s and p_e
        """
        p_s = self.convert_v(p_s)
        p_e = self.convert_v(p_e)

        if np.allclose(p_s - self.v, 0):
            if length is True:
                return True, 0, np.sqrt(np.sum(np.square(p_e-self.v)))
            else:
                return True
        elif np.allclose(p_e - self.v, 0):
            if length is True:
                return True, np.sqrt(np.sum(np.square(p_s-self.v))), 0
            else:
                return True
        else:
            l_1 = (p_e-p_s) / np.linalg.norm(p_e-p_s)
            l_2 = (p_s-p_e) / np.linalg.norm(p_s-p_e)
            c_1 = (p_e-self.v) / np.linalg.norm(p_e-self.v)
            c_2 = (p_s-self.v) / np.linalg.norm(p_s-self.v)
            if np.allclose(l_1, c_1):
                if np.allclose(l_2, c_2):
                    if length is True:
                        return (True, np.sqrt(np.sum(np.square(p_s-self.v))),
                                np.sqrt(np.sum(np.square(p_e-self.v))))
                    else:
                        return True
                return False
            else:
                return False

    def in_rect(self, p_1, p_2, p_3, p_4):
        """
        Checks if the point is inside or on one of the sides of the rectangle
        made out of the four points supplied.
        Also checks if the four points supplied actually form a rectangle.
        https://stackoverflow.com/questions/2752725/
        finding-whether-a-point-lies-inside-a-rectangle-or-not
        """
        p_1 = self.convert_v(p_1)
        p_2 = self.convert_v(p_2)
        p_3 = self.convert_v(p_3)
        p_4 = self.convert_v(p_4)

        if self.is_rect(p_1, p_2, p_3, p_4):

            if any(np.array_equal(self.v, k) for k in [p_1, p_2, p_3, p_4]):
                return True

            s_12 = p_2 - p_1
            s_23 = p_3 - p_2
            s_p1 = self.v - p_1
            s_p2 = self.v - p_2

            if all((0 <= np.dot(s_p1, s_12) <= np.dot(s_12, s_12),
                    0 <= np.dot(s_p2, s_23) <= np.dot(s_23, s_23))):
                return True
            else:
                return False
        else:
            print("check_point = %s\n"
                  "A = %s\n"
                  "B = %s\n"
                  "C = %s\n"
                  "D = %s\n" % (self, p_1, p_2, p_3, p_4))
            return False

    def inside_rect(self, p_1, p_2, p_3, p_4):
        """
        Checks if the point is inside of the rectangle
        made out of the four points supplied.
        Also checks if the four points supplied actually form a rectangle.
        https://stackoverflow.com/questions/2752725/
        finding-whether-a-point-lies-inside-a-rectangle-or-not
        """

        p_1 = self.convert_v(p_1)
        p_2 = self.convert_v(p_2)
        p_3 = self.convert_v(p_3)
        p_4 = self.convert_v(p_4)

        if self.is_rect(p_1, p_2, p_3, p_4):

            s_12 = p_2 - p_1
            s_23 = p_3 - p_2
            s_p1 = self.v - p_1
            s_p2 = self.v - p_2

            if all((0 <= np.dot(s_p1, s_12) <= np.dot(s_12, s_12),
                    0 <= np.dot(s_p2, s_23) <= np.dot(s_23, s_23))):
                if not self.on_rect(p_1, p_2, p_3, p_4):
                    return True
                else:
                    return False
            else:
                return False
        else:
            print("check_point = %s\n"
                  "A = %s\n"
                  "B = %s\n"
                  "C = %s\n"
                  "D = %s\n" % (self, p_1, p_2, p_3, p_4))
            return False

    def is_rect(self, p_1, p_2, p_3, p_4):
        """
        Checks if four points supplied form a rectangle (or a square).
        """
        p_1 = self.convert_v(p_1)
        p_2 = self.convert_v(p_2)
        p_3 = self.convert_v(p_3)
        p_4 = self.convert_v(p_4)

        p_l = [p_1, p_2, p_3, p_4]

        for pair in combinations(p_l, 2):
            if np.array_equal(pair[0], pair[1]):
                return False

        d_lst = []
        d_lst.append(np.sqrt(np.sum(np.square(p_1 - p_2))))
        d_lst.append(np.sqrt(np.sum(np.square(p_2 - p_3))))
        d_lst.append(np.sqrt(np.sum(np.square(p_3 - p_4))))
        d_lst.append(np.sqrt(np.sum(np.square(p_4 - p_1))))
        d_lst.append(np.sqrt(np.sum(np.square(p_4 - p_2))))
        d_lst.append(np.sqrt(np.sum(np.square(p_3 - p_1))))

        z = set(d_lst)

        if len(z) == 3:
            return True
        elif np.isclose(np.min(d_lst) * np.sqrt(2), np.max(d_lst)):
            return True
        else:
            return False

    def on_rect(self, p_1, p_2, p_3, p_4):
        """
        Checks if this point is on one of the sides of the rectangle formed by
        the four points supplied. If it is, it also returns the points
        comprising the side the point is on.
        Also checks if the four points supplied form a rectangle.
        """
        if self.is_rect(p_1, p_2, p_3, p_4):
            if self.on_line(p_1, p_2):
                return True, p_1, p_2
            elif self.on_line(p_2, p_3):
                return True, p_2, p_3
            elif self.on_line(p_3, p_4):
                return True, p_3, p_4
            elif self.on_line(p_4, p_1):
                return True, p_4, p_1
            else:
                return False
        else:
            False

    def in_between(self, other):
        """
        Returns a point in between this point and the other point or location
        given.
        """
        p_e = self.convert_v(other)

        p_s = self.v

        p_n = Point(*(p_s + ((p_e-p_s)/2)))

        if p_n.on_line(p_s, p_e):
            return p_n

    def intersect(self, p_1, p_2, p_3, p_4):
        """
        Determines if the lines comprised of points [p_1, p_2] and [p_3, p_4]
        intersect and returns the point where they intersect if they do.
        """
        d_1 = p_1.v - p_2.v
        d_2 = p_4.v - p_3.v
        d_3 = p_3.v - p_1.v
        if not np.allclose(np.cross(d_1, d_3), 0):
            if not np.allclose(np.cross(d_1, d_2), 0):
                if np.allclose(np.cross(d_1, d_2), 0):
                    k = np.nan_to_num(np.cross(d_1, d_3) / np.cross(d_1, d_2))
                    s = p_1.v + k * d_1
                else:
                    k = np.nan_to_num(np.cross(d_1, d_3) / np.cross(d_1, d_2))
                    s = p_1.v - k * d_1
                return True, s
            else:
                return False
        else:
            return False

    def parallel(self, p_1, p_2, p_3, p_4):
        """
        Determines if the lines comprised of points [p_1, p_2] and [p_3, p_4]
        are parallel.
        """
        d_1 = p_1.v - p_2.v
        d_2 = p_4.v - p_3.v

        if np.allclose(np.cross(d_1, d_2), 0):
            return True
        else:
            return False

    def convert_v(self, point):
        """
        returns a numpy array of the point, regardless if it is supplied as a
        point class, list or numpy array.
        """
        try:
            return point.v
        except AttributeError:
            return np.array(point)

    def distance(self, point):
        """
        Calculates the distance from the point supplied to this point.
        """
        point = self.convert_v(point)
        return np.sqrt(np.sum(np.square(self.v - point)))
