import numpy as np
import time
from edsFEM.points import Point
from edsFEM.systems import System
from edsFEM.renderers import Renderer
from edsFEM.profiles import HEB_1000, HEB_800
from edsFEM.sections import WeldedIBeam, Stiffeners

if __name__ == '__main__':

    """
    In this example a preliminary design is made for the main beams in an
    industrial floor system. The beams are loaded by big point loads from
    heavy machinery and high distributed loads are applied to the floor area.

    When the stresses in the beams are checked, it is clear that the design in
    its current state is not sufficient.
    """

    try:

        render = 1  # Render the result
        sw = 1.35  # selfweight factor

        f_DL = 1.35  # Dead-load factor
        f_LL = 1.5  # Live-load factor

        t_t = time.perf_counter()  # Timing

        #  Dimensions

        x_l = 18800  # mm

        x_1 = 9400  # mm
        x_2 = 9400  # mm

        y_l = 21200  # mm
        y_1 = 11050  # mm
        y_2 = 10150  # mm

        #  Initialize the system

        t_s = time.perf_counter()

        fem_sys = System()

        print()
        print("Initializing the system took: " +
              "{:.2f}".format(time.perf_counter()-t_s) + " s")
        print()

        t_s = time.perf_counter()

        # Beam segments

        s = 13

        # Releases

        r_1_2 = np.array([0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1])
        r_1 = np.array([0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0])
        r_2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1])
        r_0 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        #  Cyclone beams

        c_Fz = 4000e+3  # N vertical load
        c_Fh = 0.5 * c_Fz  # N Horizontal load

        c_m_prof = WeldedIBeam(HEB_1000.E, HEB_1000.G, 1800, 850, 95, 200, 30)
        c_s_prof = HEB_1000
        c_s_stiff = Stiffeners('hor', 20)
        c_s_prof.add_stiffeners(c_s_stiff)

        c_points = []

        c_points.append(Point(0, y_1, 0))
        c_points.append(Point(9270, y_1, 0))
        c_points.append(Point(17190, y_1, 0))
        c_points.append(Point(x_l, y_1, 0))

        c_points.append(Point(0, y_1-7800, 0))
        c_points.append(Point(9270, y_1-7800, 0))
        c_points.append(Point(17190, y_1-7800, 0))
        c_points.append(Point(x_l, y_1-7800, 0))

        c_beams = []
        c_loads = []

        c_beams.append(fem_sys.add_beam(c_points[0],
                                        c_points[1],
                                        profile=c_m_prof,
                                        release=r_1,
                                        segments=s))

        c_beams.append(fem_sys.add_beam(c_points[1],
                                        c_points[2],
                                        profile=c_m_prof,
                                        segments=s))

        c_loads.append(fem_sys.add_p_load(c_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * c_Fz/4,
                                          Fy=f_DL * c_Fh/4))

        c_beams.append(fem_sys.add_beam(c_points[2],
                                        c_points[3],
                                        profile=c_m_prof,
                                        release=r_2,
                                        segments=s))

        c_beams.append(fem_sys.add_beam(c_points[4],
                                        c_points[5],
                                        profile=c_m_prof,
                                        release=r_1,
                                        segments=s))

        c_beams.append(fem_sys.add_beam(c_points[5],
                                        c_points[6],
                                        profile=c_m_prof,
                                        segments=s))

        c_loads.append(fem_sys.add_p_load(c_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * c_Fz/4,
                                          Fy=-f_DL * c_Fh/4))

        c_beams.append(fem_sys.add_beam(c_points[6],
                                        c_points[7],
                                        profile=c_m_prof,
                                        release=r_2,
                                        segments=s))

        c_beams.append(fem_sys.add_beam(c_points[1],
                                        c_points[5],
                                        profile=c_s_prof,
                                        release=r_0,
                                        segments=s))

        c_loads.append(fem_sys.add_p_load(c_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * c_Fz/4,
                                          Fx=-f_DL * c_Fh/4))

        c_beams.append(fem_sys.add_beam(c_points[2],
                                        c_points[6],
                                        profile=c_s_prof,
                                        release=r_0,
                                        segments=s))

        c_loads.append(fem_sys.add_p_load(c_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * c_Fz/4,
                                          Fx=f_DL * c_Fh/4))

        #  Pyroclon beams

        p_Fz = 1200e+3  # N
        p_Fh = 0.5 * c_Fz

        p_m_prof = WeldedIBeam(HEB_1000.E, HEB_1000.G, 1800, 850, 80, 140, 30)
        p_s_prof = HEB_800

        p_points = []

        p_points.append(Point(0, 14770, 0))
        p_points.append(Point(6820, 14770, 0))
        p_points.append(Point(11700, 14770, 0))
        p_points.append(Point(12440, 14770, 0))
        p_points.append(Point(17300, 14770, 0))
        p_points.append(Point(x_l, 14770, 0))

        p_points.append(Point(0, 19860, 0))
        p_points.append(Point(6820, 19860, 0))
        p_points.append(Point(11700, 19860, 0))
        p_points.append(Point(12440, 19860, 0))
        p_points.append(Point(17300, 19860, 0))
        p_points.append(Point(x_l, 19860, 0))

        p_beams = []
        p_loads = []

        p_beams.append(fem_sys.add_beam(p_points[0],
                                        p_points[1],
                                        profile=p_m_prof,
                                        release=r_1,
                                        segments=s))

        p_beams.append(fem_sys.add_beam(p_points[1],
                                        p_points[2],
                                        profile=p_m_prof,
                                        segments=s))

        p_loads.append(fem_sys.add_p_load(p_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * p_Fz/8,
                                          Fy=-f_DL * p_Fh/8))

        p_beams.append(fem_sys.add_beam(p_points[2],
                                        p_points[3],
                                        profile=p_m_prof,
                                        segments=s))

        p_beams.append(fem_sys.add_beam(p_points[3],
                                        p_points[4],
                                        profile=p_m_prof,
                                        segments=s))

        p_loads.append(fem_sys.add_p_load(p_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * p_Fz/8,
                                          Fy=-f_DL * p_Fh/8))

        p_beams.append(fem_sys.add_beam(p_points[4],
                                        p_points[5],
                                        profile=p_m_prof,
                                        release=r_2,
                                        segments=s))

        p_beams.append(fem_sys.add_beam(p_points[6],
                                        p_points[7],
                                        profile=p_m_prof,
                                        release=r_1,
                                        segments=s))

        p_beams.append(fem_sys.add_beam(p_points[7],
                                        p_points[8],
                                        profile=p_m_prof,
                                        segments=s))

        p_loads.append(fem_sys.add_p_load(p_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * p_Fz/8,
                                          Fy=f_DL * p_Fh/8))

        p_beams.append(fem_sys.add_beam(p_points[8],
                                        p_points[9],
                                        profile=p_m_prof,
                                        segments=s))

        p_beams.append(fem_sys.add_beam(p_points[9],
                                        p_points[10],
                                        profile=p_m_prof,
                                        segments=s))

        p_loads.append(fem_sys.add_p_load(p_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * p_Fz/8,
                                          Fy=f_DL * p_Fh/8))

        p_beams.append(fem_sys.add_beam(p_points[10],
                                        p_points[11],
                                        profile=p_m_prof,
                                        release=r_2,
                                        segments=s))

        p_beams.append(fem_sys.add_beam(p_points[1],
                                        p_points[7],
                                        profile=p_s_prof,
                                        release=r_0,
                                        segments=s))

        p_loads.append(fem_sys.add_p_load(p_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * p_Fz/8,
                                          Fx=-f_DL * p_Fh/8))

        p_beams.append(fem_sys.add_beam(p_points[2],
                                        p_points[8],
                                        profile=p_s_prof,
                                        release=r_0,
                                        segments=s))

        p_loads.append(fem_sys.add_p_load(p_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * p_Fz/8,
                                          Fx=f_DL * p_Fh/8))

        p_beams.append(fem_sys.add_beam(p_points[3],
                                        p_points[9],
                                        profile=p_s_prof,
                                        release=r_0,
                                        segments=s))

        p_loads.append(fem_sys.add_p_load(p_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * p_Fz/8,
                                          Fx=-f_DL * p_Fh/8))

        p_beams.append(fem_sys.add_beam(p_points[4],
                                        p_points[10],
                                        profile=p_s_prof,
                                        release=r_0,
                                        segments=s))

        p_loads.append(fem_sys.add_p_load(p_beams[-1].nodes[int(s/2)],
                                          Fz=-f_DL * p_Fz/8,
                                          Fx=f_DL * p_Fh/8))

        #  Outside beams

        o_m_prof = WeldedIBeam(HEB_1000.E, HEB_1000.G, 1500, 400, 40, 75, 20)
        o_s_prof = HEB_800

        o_points = []

        o_points.append(Point(0, 0, 0))

        o_points.append(Point(0, y_l, 0))

        o_points.append(Point(x_1, y_l, 0))
        o_points.append(Point(x_l, y_l, 0))

        o_points.append(Point(x_l, 0, 0))
        o_points.append(Point(x_1, 0, 0))

        o_beams = []

        o_beams.append(fem_sys.add_beam(o_points[0],
                                        c_points[4],
                                        profile=o_m_prof,
                                        release=r_0,
                                        segments=s))

        fem_sys.add_support(o_beams[-1].nodes[0],
                            t_x=0, t_y=0, t_z=0,
                            r_x=False, r_y=False, r_z=False)

        o_beams.append(fem_sys.add_beam(c_points[4],
                                        c_points[0],
                                        profile=o_m_prof,
                                        release=r_0,
                                        segments=s))

        fem_sys.add_support(o_beams[-1].nodes[-1],
                            t_x=0, t_y=0, t_z=0,
                            r_x=False, r_y=False, r_z=False)

        o_beams.append(fem_sys.add_beam(c_points[0],
                                        p_points[0],
                                        profile=o_m_prof,
                                        release=r_0,
                                        segments=s))

        o_beams.append(fem_sys.add_beam(p_points[0],
                                        p_points[6],
                                        profile=o_m_prof,
                                        segments=s))

        o_beams.append(fem_sys.add_beam(p_points[6],
                                        o_points[1],
                                        profile=o_m_prof,
                                        release=r_0,
                                        segments=s))

        fem_sys.add_support(o_beams[-1].nodes[-1],
                            t_x=0, t_y=0, t_z=0,
                            r_x=False, r_y=False, r_z=False)

        o_beams.append(fem_sys.add_beam(o_points[1],
                                        o_points[2],
                                        profile=o_s_prof,
                                        release=r_0,
                                        segments=s))

        fem_sys.add_support(o_beams[-1].nodes[-1],
                            t_x=0, t_y=0, t_z=0,
                            r_x=False, r_y=False, r_z=False)

        o_beams.append(fem_sys.add_beam(o_points[2],
                                        o_points[3],
                                        profile=o_s_prof,
                                        release=r_0,
                                        segments=s))

        fem_sys.add_support(o_beams[-1].nodes[-1],
                            t_x=0, t_y=0, t_z=0,
                            r_x=False, r_y=False, r_z=False)

        o_beams.append(fem_sys.add_beam(o_points[3],
                                        p_points[11],
                                        profile=o_m_prof,
                                        release=r_0,
                                        segments=s))

        o_beams.append(fem_sys.add_beam(p_points[5],
                                        p_points[11],
                                        profile=o_m_prof,
                                        segments=s))

        o_beams.append(fem_sys.add_beam(p_points[5],
                                        c_points[3],
                                        profile=o_m_prof,
                                        release=r_0,
                                        segments=s))

        fem_sys.add_support(o_beams[-1].nodes[-1],
                            t_x=0, t_y=0, t_z=0,
                            r_x=False, r_y=False, r_z=False)

        o_beams.append(fem_sys.add_beam(c_points[3],
                                        c_points[7],
                                        profile=o_m_prof,
                                        segments=s))

        o_beams.append(fem_sys.add_beam(c_points[7],
                                        o_points[4],
                                        profile=o_m_prof,
                                        release=r_0,
                                        segments=s))

        fem_sys.add_support(o_beams[-1].nodes[-1],
                            t_x=0, t_y=0, t_z=0,
                            r_x=False, r_y=False, r_z=False)

        o_beams.append(fem_sys.add_beam(o_points[4],
                                        o_points[5],
                                        profile=o_s_prof,
                                        release=r_0,
                                        segments=s))

        fem_sys.add_support(o_beams[-1].nodes[-1],
                            t_x=0, t_y=0, t_z=0,
                            r_x=False, r_y=False, r_z=False)

        o_beams.append(fem_sys.add_beam(o_points[5],
                                        o_points[0],
                                        profile=o_s_prof,
                                        release=r_0,
                                        segments=s))

        print()
        print("Setting up took: " +
              "{:.2f}".format(time.perf_counter()-t_s) + " s")
        print()

        t_l = time.perf_counter()

        #  Add square loads

        q_live = f_LL * 5e-3  # N / mm**2
        q_dead = f_DL * (55 * 9.81 * 10e-6)  # N / mm**2 (cheq. plates)

        q = q_live + q_dead

        fem_sys.add_square_load(q, o_points[0], c_points[-1],
                                direction=[0, 0, -1],
                                auto=True)

        fem_sys.add_square_load(q, c_points[4], c_points[1],
                                direction=[0, 0, -1],
                                auto=True)

        fem_sys.add_square_load(q, c_points[-2], c_points[3],
                                direction=[0, 0, -1],
                                auto=True)

        fem_sys.add_square_load(q, c_points[0], p_points[5],
                                direction=[0, 0, -1],
                                auto=True)

        fem_sys.add_square_load(q, p_points[0], p_points[7],
                                direction=[0, 0, -1],
                                auto=True)

        fem_sys.add_square_load(q, p_points[2], p_points[9],
                                direction=[0, 0, -1],
                                auto=True)

        fem_sys.add_square_load(q, p_points[4], p_points[-1],
                                direction=[0, 0, -1],
                                auto=True)

        fem_sys.add_square_load(q, p_points[6], o_points[3],
                                direction=[0, 0, -1],
                                auto=True)

        print()
        print("Adding loads took: " +
              "{:.2f}".format(time.perf_counter()-t_l) + ' s')
        print()

        t_3 = time.perf_counter()
        fem_sys.assemble(selfweight=sw, t_test=True)
        t_4 = time.perf_counter()

        print()
        print("Total assembly took: " +
              "{:.2f}".format(t_4-t_3) + " s")
        print()

        t_s = time.perf_counter()
        fem_sys.solve(t_test=True)
        t_s = time.perf_counter() - t_s

        print()
        print("Solving took: " +
              "{:.2f}".format(t_s) + " s")
        print()

        t_p = time.perf_counter()
        fem_sys.postprocess()
        t_p = time.perf_counter() - t_p
        print()
        print("Postprocessing took: " +
              "{:.2f}".format(t_p) + " s")
        print()

        max_ym = 0
        max_zm = 0
        G_low = 0
        G_mid = 0
        G_high = 0
        G_stiff = 0

        for beam in fem_sys.beams:
            if np.max(np.abs(beam.My[:, 3])) > max_ym:
                max_ym = np.max(np.abs(beam.My[:, 3]))
            if np.max(np.abs(beam.Mz[:, 3])) > max_zm:
                max_zm = np.max(np.abs(beam.My[:, 3]))
            g_b = (beam.profile.g-beam.profile.g_stiff)*1000
            if g_b <= 40:
                G_low += (beam.profile.g-beam.profile.g_stiff) * beam.L
            elif 40 < g_b <= 100:
                G_mid += (beam.profile.g-beam.profile.g_stiff) * beam.L
            elif 100 < g_b:
                G_high += (beam.profile.g-beam.profile.g_stiff) * beam.L
            G_stiff += beam.profile.g_stiff * beam.L
        G = [G_low, G_mid, G_high, G_stiff]

        print()
        print("Beams < 40 kg/m :              " +
              "{:.2f}".format(G_low * 10**-3) + " tonnes")
        print("40 kg/m < Beams < 100 kg/m:    " +
              "{:.2f}".format(G_mid * 10**-3) + " tonnes")
        print("100 kg/m < Beams:              " +
              "{:.2f}".format(G_high * 10**-3) + " tonnes")
        print("Stiffener weight:              " +
              "{:.2f}".format(G_stiff * 10**-3) + " tonnes")

        print("Total weight:                  " +
              "{:.2f}".format(np.sum(G) * 10**-3) + " tonnes")
        print()

        if render == 1:

            t_r = time.perf_counter()

            renderbox = Renderer()

            renderbox.render_origin(scale='mm')

            if fem_sys.structurenodes is not None:
                for node in fem_sys.structurenodes:
                    renderbox.render_node(node, scale='mm')
            else:
                if fem_sys.nodes is not None:
                    for node in fem_sys.nodes:
                        renderbox.render_node(node, scale='mm')

            if fem_sys.beams is not None:
                for beam in fem_sys.beams:
                    renderbox.render_beam(beam, scale='mm')

            renderbox.make_button('supports', [fem_sys.nodes, 'mm', 8])
            renderbox.make_button('beams', [fem_sys.beams, 'mm', None])
            renderbox.make_button('releases', [fem_sys.beams, 'mm', 6])
            renderbox.make_button('sigma', [fem_sys.beams, 'mm', 235])
            renderbox.make_button('shear', [fem_sys.beams, 'mm', 235])
            renderbox.make_button('v_mises', [fem_sys.beams, 'mm', 235])
            renderbox.make_button('displacements', [fem_sys.beams, 'mm', 6])
            renderbox.make_button('profiles', [fem_sys.beams, 'mm'])
            renderbox.make_button('moments', [fem_sys.beams, 'mm', 1e-6])
            renderbox.make_button('n_forces', [fem_sys.beams, 'mm', 1e-3])
            renderbox.make_button('v_forces', [fem_sys.beams, 'mm', 1e-3])
            renderbox.make_button('p_loads', [fem_sys.beams, 'mm', 100e+3])
            renderbox.make_button('m_loads', [fem_sys.beams, 'mm', 5e+4])

            print()
            print("Rendering took: " +
                  "{:.2f}".format(time.perf_counter() - t_r) + " s")
            print()

            print()
            print("Total took: " +
                  "{:.2f}".format(time.perf_counter() - t_t) + " s")
            print()

            renderbox.run()

        else:

            print()
            print("Total took: " +
                  "{:.2f}".format(time.perf_counter() - t_t) + " s")
            print()

    except:
        import sys
        print(sys.exc_info()[0])
        import traceback
        print(traceback.format_exc())
    finally:
        input("Press enter to exit")
