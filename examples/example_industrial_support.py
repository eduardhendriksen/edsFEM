import copy
import numpy as np
from edsFEM.points import Point
from edsFEM.systems import System
from edsFEM.renderers import Renderer
from edsFEM.profiles import IPE, HE
from edsFEM.sections import Stiffeners

if __name__ == '__main__':

    """
    In this example a design is made for a supporting construction in an
    industrial setting. The maximal vertical deflection for the top of the
    support construction is set to be at most 5mm.
    The designer has chosen to use the same profile for all beams in this
    example.
    """

    Factor = 1.35

    F_8 = 360e+3  # N  Vertical load
    F_9 = 610e+3  # N  Vertical load

    n_8_abst = 2  # number of support plates
    n_9_abst = 4  # number of support plates

    stiff_9 = 15  # mm Stiffener thickness
    stiff_8 = 7  # mm Stiffener thickness

    s = 12  # No. of beam segments

    prof_same = 'HEB 600'  # Profile used

    render = 1  # Toggles rendering

    try:

        # Geometry

        p_0 = Point(0, 0, 4600)
        p_1 = Point(1110, 0, 4600)
        p_2 = Point(4410, 1000, 4600)
        p_3 = Point(3460, 5125, 4600)
        p_4 = Point(1110, 5125, 4600)
        p_5 = Point(-956, 4500, 4600)

        p_6 = Point(4600, 0, 4600)
        p_7 = Point(5935, 0, 4600)
        p_8 = Point(5935, 5125, 4600)

        p_9 = Point(-1110, 5125, 4600)

        Fz_9 = (Factor * F_9) / n_9_abst
        Fh_9 = Fz_9 * 0.5

        #  Profiles

        stiff_9 = [Stiffeners('hor', stiff_9)]  # Stiffeners applied

        p_s_9 = copy.deepcopy([p for p in IPE+HE if
                               p.designation == prof_same][0])
        p_s_9.add_stiffeners(stiff_9)
        prof_9 = [p for p in IPE+HE if p.designation == prof_same][0]

        #  Initiate system

        fem_sys = System()

        #  Add beams

        m_b_9 = []

        m_b_9.append(fem_sys.add_beam(p_0, p_1, profile=p_s_9))
        m_b_9.append(fem_sys.add_beam(p_1, p_2, profile=p_s_9))
        m_b_9.append(fem_sys.add_beam(p_2, p_3, profile=p_s_9))
        m_b_9.append(fem_sys.add_beam(p_3, p_4, profile=p_s_9))
        m_b_9.append(fem_sys.add_beam(p_4, p_5, profile=p_s_9))
        m_b_9.append(fem_sys.add_beam(p_5, p_0, profile=p_s_9))

        fem_sys.add_beam(p_1, p_6, profile=prof_9)
        fem_sys.add_beam(p_6, p_2, profile=prof_9)
        fem_sys.add_beam(p_6, p_7, profile=prof_9)
        fem_sys.add_beam(p_7, p_8, profile=prof_9)
        fem_sys.add_beam(p_8, p_3, profile=prof_9)
        fem_sys.add_beam(p_4, p_9, profile=prof_9)
        fem_sys.add_beam(p_9, p_5, profile=prof_9)

        fem_sys.add_p_load(fem_sys.beams[1].nodes[4], Fz=-Fz_9, Fy=-Fh_9)
        fem_sys.add_p_load(fem_sys.beams[2].nodes[6], Fz=-Fz_9, Fx=Fh_9)
        fem_sys.add_p_load(fem_sys.beams[4].nodes[0], Fz=-Fz_9, Fy=Fh_9)
        fem_sys.add_p_load(fem_sys.beams[5].nodes[int(s/2)],
                           Fz=-Fz_9, Fx=-Fh_9)

        p_10 = Point(1110, 0, 0)
        p_11 = Point(5935, 0, 0)
        p_12 = Point(5935, 5125, 0)
        p_13 = Point(1110, 5125, 0)

        Fz_8 = (Factor * F_8) / n_8_abst
        Fh_8 = Fz_8 * 0.5

        #  Profiles

        stiff_8 = [Stiffeners('hor', stiff_8)]

        p_s_8 = copy.deepcopy([p for p in IPE+HE if
                               p.designation == prof_same][0])
        p_s_8.add_stiffeners(stiff_8)

        prof_8 = [p for p in IPE+HE if p.designation == prof_same][0]

        #  Add beams

        m_b_8 = []

        m_b_8.append(fem_sys.add_beam(p_10, p_11, profile=p_s_8))

        fem_sys.add_p_load(fem_sys.beams[-1].nodes[int(s/2)],
                           Fz=-Fz_8, Fy=-Fh_8)

        m_b_8.append(fem_sys.add_beam(p_12, p_13, profile=p_s_8))

        fem_sys.add_p_load(fem_sys.beams[-1].nodes[int(s/2)],
                           Fz=-Fz_8, Fy=Fh_8)

        fem_sys.add_beam(p_11, p_12, profile=prof_8)
        fem_sys.add_beam(p_13, p_10, profile=prof_8)

        fem_sys.add_beam(p_1, p_10, profile=prof_8)
        fem_sys.add_beam(p_7, p_11, profile=prof_8)
        fem_sys.add_beam(p_8, p_12, profile=prof_8)
        fem_sys.add_beam(p_4, p_13, profile=prof_8)

        p_14 = Point(500, 0, 0)
        p_15 = Point(8000, 0, 0)
        p_16 = Point(8000, 5125, 0)
        p_17 = Point(500, 5125, 0)

        fem_sys.add_beam(p_14, p_10, profile=prof_8)
        fem_sys.add_beam(p_15, p_11, profile=prof_8)
        fem_sys.add_beam(p_16, p_12, profile=prof_8)
        fem_sys.add_beam(p_17, p_13, profile=prof_8)

        fem_sys.add_support(fem_sys.beams[-4].nodes[0], t_x=0, t_y=0, t_z=0,
                            r_x=None, r_y=None, r_z=None)
        fem_sys.add_support(fem_sys.beams[-3].nodes[0], t_x=0, t_y=0, t_z=0,
                            r_x=None, r_y=None, r_z=None)
        fem_sys.add_support(fem_sys.beams[-2].nodes[0], t_x=0, t_y=0, t_z=0,
                            r_x=None, r_y=None, r_z=None)
        fem_sys.add_support(fem_sys.beams[-1].nodes[0], t_x=0, t_y=0, t_z=0,
                            r_x=None, r_y=None, r_z=None)

        fem_sys.assemble(selfweight=Factor)
        fem_sys.solve()
        fem_sys.postprocess()

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
            renderbox.make_button('displacements', [fem_sys.beams, 'mm', 12])
            renderbox.make_button('profiles', [fem_sys.beams, 'mm'])
            renderbox.make_button('moments', [fem_sys.beams, 'mm', 2e-6])
            renderbox.make_button('n_forces', [fem_sys.beams, 'mm', 2e-3])
            renderbox.make_button('v_forces', [fem_sys.beams, 'mm', 1e-3])
            renderbox.make_button('p_loads', [fem_sys.beams, 'mm', 250e+2])
            renderbox.make_button('m_loads', [fem_sys.beams, 'mm', 5e+4])

            renderbox.run()

    except:
        import sys
        print(sys.exc_info()[0])
        import traceback
        print(traceback.format_exc())
    finally:
        input("Press enter to exit")
