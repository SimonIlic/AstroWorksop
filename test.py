import rebound
import numpy as np

sim = rebound.Simulation.from_file("solar_system.bin")

fig = rebound.OrbitPlot(sim)
fig.savefig("image.png") # save figure to file
fig.show() # show fipigure on screen
print("done")