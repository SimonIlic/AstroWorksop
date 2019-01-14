import rebound
import matplotlib.pyplot as plt

sim = rebound.Simulation()
sim.add(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn",
         "Uranus", "Neptune"])

fig = rebound.OrbitPlot(sim)
fig.savefig("circumbinary_j.png")
plt.show()
