# Jasper Lankhorst
# Program to make a gif of a simulation of planets

import imageio
import matplotlib.pyplot as plt
import numpy as np
import rebound

# start simulation of solar system
sim = rebound.Simulation.from_file("solar_system.bin")
rebound.OrbitPlot(sim)

# add planet to fly by and integrate
sim.add(m=1,x=40, y=-1000, vy=1.2)
sim.integrate(sim.t+600)

# center plot at sun
sim.move_to_hel()

# make r plots
r = 3
images = np.array([])
for i in range(r):
    sim.integrate(sim.t+10)
    fig = rebound.OrbitPlot(sim,color=True,unitlabel="[AU]", lim=250.,
        show_orbit=True)
    print ("%s of 300." %i)
    plt.savefig("plotje_%s.png" %i)
    np.append(images, "plotje_%s.png" %i)
    plt.close('all')
imageio.mimsave('fly_by_gifje.gif', images)
