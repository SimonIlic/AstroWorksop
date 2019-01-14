import rebound
import matplotlib.pyplot as plt

sim = rebound.Simulation()
sim.add(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn",
         "Uranus", "Neptune"], date="99999-01-01 06:25")

fig = rebound.OrbitPlot(sim)
fig.savefig("circumbinary_j.png")
plt.show()
