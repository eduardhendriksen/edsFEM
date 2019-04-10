import time
import numpy as np
from edsFEM.points import Point
from edsFEM.systems import System
from edsFEM.renderers import Renderer
from edsFEM.profiles import HEA_450
from edsFEM.optimizers import BeamOptimizer

if __name__ == '__main__':

    """
    In this example a simple floor system is loaded automatically using the
    square load distribution function to load the beams. Subsequently the beams
    are roughly 'optimized' based on a heavily simplified Von Mises stress
    calculation in multiple loops.
    """

    try:

        render = 1  # Toggles rendering
        factor_LL = 1.5  # Live load factor
        factor_DL = 1.35  # Dead load factor

        sw = factor_DL  # self-weight factor

        s = 13  # No. of beam segments

        # End releases

        end_releases = np.array([0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1])
        r_1 = np.array([0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0])
        r_2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1])

        t_t = time.perf_counter()  # timing

        #  Dimensions

        l_1 = 11500  # mm
        l_2 = 10300  # mm
        h_1 = 12000  # mm
        h_2 = 20900  # mm
        h_3 = 26100  # mm

        p_0 = Point(0, 0, 0)
        p_1 = Point(l_1/2, 0, 0)
        p_2 = Point(l_1, 0, 0)
        p_3 = Point(0, l_2, 0)
        p_4 = Point(l_1/2, l_2/2, 0)
        p_5 = Point(3*l_1/4, l_2/2, 0)
        p_6 = Point(l_1/2, l_2, 0)
        p_7 = Point(3*l_1/4, l_2, 0)
        p_8 = Point(l_1, l_2/2, 0)
        p_9 = Point(l_1, l_2, 0)

        #  Initialize the system

        t_1 = time.perf_counter()
        fem_sys = System()

        #  Add beams

        beams = []

        beams.append(fem_sys.add_beam(p_0, p_2, profile=HEA_450,
                                      release=end_releases, segments=s))
        fem_sys.add_support(beams[-1].nodes[0], t_x=0, t_y=0, t_z=0,
                            r_x=1, r_y=1, r_z=1)
        fem_sys.add_support(beams[-1].nodes[-1], t_x=0, t_y=0, t_z=0,
                            r_x=1, r_y=1, r_z=1)

        beams.append(fem_sys.add_beam(p_0, p_3, profile=HEA_450,
                                      release=end_releases, segments=s))
        fem_sys.add_support(beams[-1].nodes[0], t_x=0, t_y=0, t_z=0,
                            r_x=1, r_y=1, r_z=1)

        beams.append(fem_sys.add_beam(p_3, p_9, profile=HEA_450,
                                      release=end_releases, segments=s))
        fem_sys.add_support(beams[-1].nodes[0], t_x=0, t_y=0, t_z=0,
                            r_x=1, r_y=1, r_z=1)

        beams.append(fem_sys.add_beam(p_1, p_4, profile=HEA_450,
                                      release=r_1, segments=s))
        fem_sys.add_support(beams[-1].nodes[0], t_x=0, t_y=0, t_z=0,
                            r_x=1, r_y=1, r_z=1)

        beams.append(fem_sys.add_beam(p_4, p_6, profile=HEA_450,
                                      release=r_2, segments=s))

        fem_sys.add_support(beams[-1].nodes[-1], t_x=0, t_y=0, t_z=0,
                            r_x=1, r_y=1, r_z=1)

        beams.append(fem_sys.add_beam(p_4, p_5, profile=HEA_450,
                                      release=r_1, segments=s))

        beams.append(fem_sys.add_beam(p_5, p_7, profile=HEA_450,
                                      release=end_releases, segments=s))

        beams.append(fem_sys.add_beam(p_5, p_8, profile=HEA_450,
                                      release=r_2, segments=s))

        beams.append(fem_sys.add_beam(p_2, p_8, profile=HEA_450,
                                      release=r_1, segments=s))
        fem_sys.add_support(beams[-1].nodes[0], t_x=0, t_y=0, t_z=0,
                            r_x=1, r_y=1, r_z=1)

        beams.append(fem_sys.add_beam(p_8, p_9, profile=HEA_450,
                                      release=r_2, segments=s))
        fem_sys.add_support(beams[-1].nodes[-1], t_x=0, t_y=0, t_z=0,
                            r_x=1, r_y=1, r_z=1)

        print()
        print("Setting up took: " +
              "{:.2f}".format(time.perf_counter()-t_1) + " s")
        print()

        t_l = time.perf_counter()

        #  Add loads

        F_filter = 133333  # N

        q_live = 1.5 * 5e-3  # N / mm**2
        q_dead = 1.35 * (55 * 9.81 * 10e-6)  # N / mm**2 (cheq. plates)

        q = q_live + q_dead

        fem_sys.add_square_load(q, p_0, p_9,
                                direction=[0, 0, -1], name='1st floor live 1',
                                auto=True)

        print()
        print("Adding loads took: " +
              "{:.2f}".format(time.perf_counter()-t_l) + " s")
        print()

        t_a = time.perf_counter()
        fem_sys.assemble(selfweight=sw, t_test=True)

        print()
        print("Total assembly took: " +
              "{:.2f}".format(time.perf_counter()-t_a) + " s")
        print()

        t_s = time.perf_counter()
        fem_sys.solve()

        print()
        print("Solving took: " +
              "{:.2f}".format(time.perf_counter()-t_s) + " s")
        print()

        t_p = time.perf_counter()
        fem_sys.postprocess()

        print()
        print("Postprocessing took: " +
              "{:.2f}".format(time.perf_counter()-t_p) + " s")
        print()

        t_o_ = time.perf_counter()
        optimizer = BeamOptimizer(fem_sys)
        optimizer.standard_sections(s_max=235)

        print()
        print("Finding new beams took: " +
              "{:.2f}".format(time.perf_counter()-t_o_) + " s")
        print()

        t_w = time.perf_counter()
        fem_sys.clear_weight()

        print()
        print("Clearing weight took: " +
              "{:.2f}".format(time.perf_counter()-t_w) + " s")
        print()

        t_a = time.perf_counter()
        fem_sys.assemble(selfweight=sw, t_test=True)

        print()
        print("Total assembly took: " +
              "{:.2f}".format(time.perf_counter()-t_a) + " s")
        print()

        t_s = time.perf_counter()
        fem_sys.solve()

        print()
        print("Solving took: " +
              "{:.2f}".format(time.perf_counter()-t_s) + " s")
        print()

        t_p = time.perf_counter()
        fem_sys.postprocess()
        print()
        print("Postprocessing took: " +
              "{:.2f}".format(time.perf_counter()-t_p) + " s")

        t_o = time.perf_counter()
        optimizer = BeamOptimizer(fem_sys)
        optimizer.standard_sections(s_max=235)

        print()
        print("Finding new beams took: " +
              "{:.2f}".format(time.perf_counter()-t_o) + " s")
        print()

        t_w = time.perf_counter()
        fem_sys.clear_weight()

        print()
        print("Clearing weight took: " +
              "{:.2f}".format(time.perf_counter()-t_w) + " s")
        print()

        t_a = time.perf_counter()
        fem_sys.assemble(selfweight=sw, t_test=True)

        print()
        print("Total assembly took: " +
              "{:.2f}".format(time.perf_counter()-t_a) + " s")
        print()

        t_s = time.perf_counter()
        fem_sys.solve()

        print()
        print("Solving took: " +
              "{:.2f}".format(time.perf_counter()-t_s) + " s")
        print()

        t_p = time.perf_counter()
        fem_sys.postprocess()
        print()
        print("Postprocessing took: " +
              "{:.2f}".format(time.perf_counter()-t_p) + " s")
        print()

        print()
        print("Optimization loops took: " +
              "{:.2f}".format(time.perf_counter() - t_o_) + " s")
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
            renderbox.make_button('moments', [fem_sys.beams, 'mm', 0.5e-5])
            renderbox.make_button('n_forces', [fem_sys.beams, 'mm', 1e-3])
            renderbox.make_button('v_forces', [fem_sys.beams, 'mm', 1e-3])
            renderbox.make_button('p_loads', [fem_sys.beams, 'mm', 1e+3])
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
