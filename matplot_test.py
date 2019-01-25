# matplot_test.py
# Jasper Lankhorst
# Program to make a gif of a simulation of planets

import imageio
import matplotlib
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

# make r plots, append to list to create gif
r = 3
imagelist = []
for i in range(r):
    sim.integrate(sim.t+10)
    fig, ax = plt.subplots()
    ax = rebound.OrbitPlot(sim,color=True,unitlabel="[AU]", lim=250.,
        show_orbit=True)
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height() + (3,))

    print(f"{i + 1} of {r}.")
    imagelist.append(image)
    plt.close('all')

with imageio.get_writer('./fly_by_gifje.gif', mode='i') as writer:
    for image in imagelist:
        writer.append_data(image)
