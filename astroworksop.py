# astroworksop.py
# programma om stabiliteit van planetenbanen te controleren
# Jasper Lankhorst
# Workshop Astronomy

import rebound
import numpy as np

from IPython.display import display, clear_output
import matplotlib.pyplot as plt

sim = rebound.Simulation()
particle_names = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune"]
# we use the NASA horizon database to look up the Sun and planets
sim.add(particle_names)

# let's give all the particles a unique hash (based on its name)
for i, particle in enumerate(sim.particles):
    particle.hash = particle_names[i]
sim.status()
sim.save("solar_system_outer_planets.bin")
sim = rebound.Simulation.from_file("solar_system_outer_planets.bin")


def simulate_fly_by(sim, intruder, visualize=False):
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
            #print(error)
            for i, particle in enumerate(sim.particles):
                distance = np.linalg.norm(particle.xyz)
                if distance > sim.exit_max_distance:
                    #print("Removed", i, str(particle.hash))
                    sim.remove(hash=particle.hash)
                    sim.move_to_com()

            return sim


def calc_escape_velocity(sim, particle):
    #sim.move_to_hel()

    r = np.linalg.norm(particle.xyz)
    G = sim.G
    m = sim.particles[0].m

    return np.sqrt(2 * G * m / r)


def strong_regime(resolution=100, n_trials=50):
    print("Starting strong regime simulation with resolution {}, {} trials each...".format(resolution, n_trials))
    xs = np.linspace(1, 50, resolution)
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
            for particle in sim.particles:
                v = np.linalg.norm(particle.vxyz)
                v_esc = calc_escape_velocity(sim, particle)
                if v > v_esc:
                    eject_count += 1
                    break
        print("Detected", eject_count, "ejections out of", n_trials, "trials.")
        f_eject[i] = eject_count / n_trials
        print(f_eject[i])


    return (xs, f_eject)


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
#     Temporary simulation, adding sun and the particle we want the orbit from
    temp_sim = rebound.Simulation()
    temp_sim.add(simulation.particles[0])
    temp_sim.add(particle)
    particle = temp_sim.particles[1]
#     Integrating over exactly one orbit
    while total_time < period:
        temp_sim.integrate(temp_sim.t+step_size)
        total_time += step_size
        locations.append(particle.xyz)

    return np.array(locations)


def check_orbit_crossing(simulation):
    """
    Checks in a simulation whether any orbits cross.
    """

#     Creating and saving lists with points on orbits
    locationslist = []
    for i, particle in enumerate(simulation.particles[1:]):
        orbit = particle.calculate_orbit()
        step_size = orbit.P * orbit.rhill / (2 * np.pi * orbit.a)
        locationslist.append(orbit_list(simulation,
                                        orbit.P, particle, step_size))

#     creating distance matrix
    for i, loc1 in enumerate(locationslist):
        for j, loc2 in enumerate(locationslist[i+1:]):
            dist_mat = spatial.distance_matrix(loc1, loc2)
            if dist_mat[np.where(dist_mat < mutual_rhill(simulation.particles[i+1],
                                                         simulation.particles[j+i+2]))].size > 0:
                print(f"Planet {i+1} and {i+j+2} (counting from star) will collide!")
                return True

    return False


def analyze_stability(sim):

    if check_immediate_ejection(sim) == True:
        return False

    elif check_orbit_crossing(sim) == True:
        return False

    elif check_kozai(sim) == True:
        return False

    elif check_AMD(sim) == True:
        return False

    else:
        return True


if __name__ == "__main__":
    analyze_stability(sim)
