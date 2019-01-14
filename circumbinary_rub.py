import rebound
import matplotlib.pyplot as plt

sim = rebound.Simulation() #initiate simulation
sim.add(m=1.0) #Star A
sim.add(m=1.0, a=1.0) # Star B
# add circumbinary planet ABb
sim.add(a=2.0)
sim.add(a=4.0)

#integrate(200)
fig = rebound.OrbitPlot(sim)
fig.savefig("circumbinary_rub.png") # save figure to file
fig.show() # show fipigure on screen

particles = sim.particles
star_A = particles[0]
star_B = particles[1]

print(star_B.vx) #

plt.show()
