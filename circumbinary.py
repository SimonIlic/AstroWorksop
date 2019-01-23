import rebound
import matplotlib.pyplot as plt

sim = rebound.Simulation() #initiate simulation
sim.add(m=1.0) #Star A
sim.add(m=1.0, a=1.0) # Star B
# add circumbinary planet ABb
sim.add(a=2.0)

fig = rebound.OrbitPlot(sim)
fig.savefig("circumbinary.png") # save figure to file
fig.show() # show fipigure on screen
