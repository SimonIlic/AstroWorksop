# Jasper Lankhorst
# PLotje bouwen voor strong regime
# Workshop astronomy 2019
# plot_strong_regime.py

import rebound
import numpy as np
from scipy import spatial

from IPython.display import display, clear_output
import matplotlib.pyplot as plt


def calc_escape_velocity(sim, particle):
    #sim.move_to_hel()

    r = np.linalg.norm(particle.xyz)
    G = sim.G
    m = sim.particles[0].m

    return np.sqrt(2 * G * m / r)


def check_immediate_ejection(sim):
    # move to Sun frame
    sim.move_to_hel()

    # calculate velocity of each particle and compare to escpae velocity
    for particle in sim.particles[1:]:
        v = np.linalg.norm(particle.vxyz)
        v_esc = calc_escape_velocity(sim, particle)
        if v >= v_esc:
            sim.move_to_com()
            return True

    sim.move_to_com()

    return False


def randomize_sim(sim):
    sim.integrate(np.random.random()*10**3)
    return sim


def simulate_fly_by(sim, intruder, visualize=False):

    intruder.hash = "intruder"
    sim.add(intruder)

    intruder_distance = np.linalg.norm(sim.particles["intruder"].xyz)
    sim.exit_max_distance = intruder_distance*1.01

    while True:
        try:
            sim.integrate(sim.t+10)

            # test_sim = sim.copy()
            # test_sim.remove(hash="intruder")
            # if check_immediate_ejection(test_sim):
            #     sim_results["v_escapes"].append(sim.t)
            #
            # if check_orbit_crossing(test_sim):
            #     sim_results["orbit_crossing"].append(sim.t)
            #
            # if check_kozai(test_sim):
            #     sim_results["kozai"].append(sim.t)
            #
            # if visualize:
            #     fig = rebound.OrbitPlot(sim,color=True,unitlabel="[AU]")
            #     display(fig)
            #     plt.close(fig)
            #     clear_output(wait=True)

        except rebound.Escape as error:
            # if type(rebound.Escape()) == type(error):
            #     print("joe")
            # remove intruder
            sim.remove(hash="intruder")
            sim.move_to_com()

            return sim


def strong_regime(resolution, n_trials):
    print(f"Starting strong regime simulation with resolution {resolution}, "
          f"{n_trials} trials each...")
    xs = np.linspace(1, 100, resolution)
    f_eject = np.ones(resolution)

    for i, x in enumerate(xs):
        print("Running r_min =", x)
        eject_count = 0.

        # run n_trials trials detecting ejection directly after fly-by
        for j in range(n_trials):
            # get a fresh simulation
            sim = rebound.Simulation.from_file("solar_system_outer_planets.bin")
            sim = randomize_sim(sim)

            intruder = rebound.Particle(m=1.,x=x,y=-1000.,vy=2.)

            sim = simulate_fly_by(sim, intruder)

            sim.move_to_hel()
            for particle in sim.particles[1:]:
                v = np.linalg.norm(particle.vxyz)
                v_esc = calc_escape_velocity(sim, particle)
                if v > v_esc:
                    eject_count += 1
                    break
        print("Detected", eject_count, "ejections out of", n_trials, "trials.")
        f_eject[i] = eject_count / n_trials
        print(f_eject[i])


    return (xs, f_eject)


xs, f_eject = strong_regime(resolution=100, n_trials=1000)

plt.plot(xs, 1-f_eject)
plt.xlabel("Periastron (AU)")
plt.ylabel("Fraction of planets that immediately escaped")
plt.xlim([0,100])
plt.show()
