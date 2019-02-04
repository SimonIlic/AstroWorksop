# astroworksop.py
# programma om stabiliteit van planetenbanen te controleren
# Jasper Lankhorst
# Workshop Astronomy

import random
import rebound
import numpy as np
from IPython.display import display, clear_output
import matplotlib.pyplot as plt
from scipy import spatial



def simulate_fly_by(sim, intruder, visualize=False):
    """
    Simulates what happens when a star flies by a planetary system.
    """
    intruder.hash = "intruder"

    sim.add(intruder)

    intruder_distance = np.linalg.norm(sim.particles["intruder"].xyz)
    sim.exit_max_distance = intruder_distance*1.01

    while True:
        try:
            sim.integrate(sim.t+5)

            if visualize:
                fig = rebound.OrbitPlot(sim,color=True,unitlabel="[AU]")
                display(fig)
                plt.close(fig)
                clear_output(wait=True)

        except rebound.Escape as error:
            # print(error)
            for i, particle in enumerate(sim.particles):
                distance = np.linalg.norm(particle.xyz)
                if distance > sim.exit_max_distance:
                    # print("Removed", i, str(particle.hash))
                    sim.remove(hash=particle.hash)
                    sim.move_to_com()

            return sim


def calc_escape_velocity(sim, particle):
    """
    For a given simulation and planets, calculates its escape velocity.
    """

    r = np.linalg.norm(particle.xyz)
    G = sim.G
    m = sim.particles[0].m

    return np.sqrt(2 * G * m / r)


def evolve_system(sim, t):
    sim.exit_max_distance = 500

    close_encounters = []
    ejections = []
    # set planet radii to their hill spheres
    for planet in sim.particles[1:]:
        planet.r = hill_radius(planet, sim.particles[0])

    sim.collision = "direct"

    end_time = sim.t + t
    while sim.t < end_time:
        try:
            sim.integrate(end_time)
        except (rebound.Collision, rebound.Escape) as error:
            sim.status()
            if type(error) == type(rebound.Collision()):
                print(error, sim.t)
                collided = []
                for particle in sim.particles:
                    if particle.lastcollision == sim.t:
                        collided.append(particle.index)

                planet_1 = sim.particles[collided[0]]
                planet_2 = sim.particles[collided[1]]

                d = np.linalg.norm(np.array(planet_1.xyz) - np.array(planet_2.xyz))
                print(planet_1.index, planet_2.index, "close encounter. distance:", d)

                resolve_collision(sim)
                close_encounters.append(((planet_1.index, planet_2.index), d, sim.t))

            else:
                print(error)
                out_of_bounds = []
                for i, particle in enumerate(sim.particles):
                    distance = np.linalg.norm(particle.xyz)
                    if distance > sim.exit_max_distance:
                        print("Removed", particle.index, str(particle.hash))
                        out_of_bounds.append(particle.hash)
                        ejections.append((particle.index, particle.xyz, sim.t))

                for hsh in out_of_bounds:
                        sim.remove(hash=hsh)
                        sim.move_to_com()

    return (sim, close_encounters, ejections)


def strong_regime(resolution, n_trials):
    """
    Simulates a fly-by and calculates whether the system is stable or not.
    """

    print(f"Starting strong regime simulation with resolution {resolution}, "
          f"{n_trials} trials each...")
    xs = np.linspace(1, 50, resolution)
    f_eject = np.ones(resolution)
    f_stable = np.ones(resolution)

    for i, x in enumerate(xs):
        print(f"Running r_min = {x} at {round(100.*x/len(xs), 1)}%")
        eject_count = 0.
        stable_count = 0.

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

            stable = analyze_stability(sim)
            if stable:
                stable_count += 1
        print("Detected", eject_count, "ejections out of", n_trials, "trials.")
        if n_trials > eject_count:
            percentage = 100 * stable_count / (n_trials-eject_count)
            print(eject_count, n_trials, stable_count)
            print(f"Of the {n_trials - eject_count} systems left, "
                  f"{percentage}% was not long term unstable.")
        f_eject[i] = eject_count / n_trials
        f_stable[i] = stable_count / n_trials

    return (xs, f_eject, f_stable)


def mutual_rhill(p1, p2):
    """
    Calculates mutual Hill radius of particle 1 and 2.
    """
    rhill_m = (p1.a + p2.a) / 2. * ((p1.m + p2.m) / 3.)**(1/3.)
    return rhill_m


def orbit_list(simulation, period, particle, step_size):
    """
    Creates list of points on an orbit.
    """

    locations = []
    total_time = 0
    # Temporary simulation, adding sun and the particle we want the orbit from
    temp_sim = rebound.Simulation()
    temp_sim.add(simulation.particles[0])
    temp_sim.add(particle)
    particle = temp_sim.particles[1]
    # Integrating over exactly one orbit
    while total_time < period:
        temp_sim.integrate(temp_sim.t+step_size)
        total_time += step_size
        locations.append(particle.xyz)

    return np.array(locations)


def check_orbit_crossing(simulation):
    """
    Checks in a simulation whether any orbits cross.
    """

    # Creating and saving lists with points on orbits
    locationslist = []
    for i, particle in enumerate(simulation.particles[1:]):
        orbit = particle.calculate_orbit()
        step_size = orbit.P * orbit.rhill / (2 * np.pi * orbit.a)
        locationslist.append(orbit_list(simulation,
                                        orbit.P, particle, step_size))

    # creating distance matrix
    for i, loc1 in enumerate(locationslist):
        for j, loc2 in enumerate(locationslist[i+1:]):
            dist_mat = spatial.distance_matrix(loc1, loc2)
            if dist_mat[np.where(dist_mat < mutual_rhill(simulation.particles[i+1],
                    simulation.particles[j+i+2]))].size > 0:
                return True

    return False


def check_immediate_ejection(sim):
    """
    Checks whether there is a planet in the simulation with v > v_escape.
    """

    # move to Sun frame
    sim.move_to_hel()

    # calculate velocity of each particle and compare to escpae velocity
    for particle in sim.particles[1:]:
        v = np.linalg.norm(particle.vxyz)
        v_esc = calc_escape_velocity(sim, particle)
        if v >= v_esc:
            return True

    return False


def check_kozai(sim):
    """
    Checks whether the kozai mechanism is happening in a simulation.
    """

    # compare all particles except the star
    for i, particle_1 in enumerate(sim.particles[1:]):
        for j, particle_2 in enumerate(sim.particles[i+2:]):
            # calculate mutual inclination. defined as difference in inclination between two orbits
            mutual_inclination = abs(particle_1.inc - particle_2.inc)
            # check if mutual inclination is between 39.2 degrees and 140.2 degrees
            if 0.684 <  mutual_inclination and mutual_inclination < 2.46:
                return True

    return False



def randomize_sim(sim):
    """
    Integrates simulation for any number of time between 0 and 999.
    """
    sim.integrate(random.randint(0, 999))
    return sim


def analyze_stability(sim):

    if check_immediate_ejection(sim) == True:
        return False

    elif check_orbit_crossing(sim) == True:
        return False

    elif check_kozai(sim) == True:
        return False

    else:
        return True


def plot_hist_first_event(time_yrs, sim_list):
    """
    Plots some stuff from the simulation given to the function in a dictionary.
    sim_list structure: [{sim: rebound.Simulation,
    intruder: rebound.Particle,
    close_encounters: array of times,
    escapes: array of times,
    kozai: array of times,
    orbit_crossing: array of times,
    v_escapes: float}, ...]
    """

    # Dict to add all lists to to make plotting easier
    lists_dict = {}

    # Create lists of first kozai measurement for inc = 90 and inc = 0
    kozai_first_90 = []
    kozai_first_0 = []
    for i, sim_dict in enumerate(sim_list):
        if sim_dict[kozai]:
            if sim_dict[intruder].z > 0:
                kozai_first_90.append(sim_dict[kozai][0])
            if sim_dict[intruder].z == 0:
                kozai_first_0.append(sim_dict[kozai][0])

    lists_dict["kozai_90"] = np.array(kozai_first_90)
    lists_dict["kozai_0"] = np.array(kozai_first_0)

    # Same for orbit crossing
    orbit_first_90 = []
    orbit_first_0 = []
    for i, sim_dict in enumerate(sim_list):
        if sim_dict[orbit_crossing]:
            if sim_dict[intruder].z > 0:
                orbit_first_90.append(sim_dict[orbit_crossing][0])
            if sim_dict[intruder].z == 0:
                orbit_first_0.append(sim_dict[orbit_crossing][0])

    lists_dict["orbit_90"] = orbit_first_90
    lists_dict["orbit_0"] = orbit_first_0


    # Plot that hist
    for i, key in enumerate(lists_dict.keys()):
        plt.hist(lists_dict[key], bins=100)
        plt.ylabel("Frequency")
        plt.xlabel("Time (yr/2pi)")
        plt.savefig(f"plot, {list_dict[key]}, histogram")




if __name__ == "__main__":
