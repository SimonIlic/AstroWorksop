import rebound
import matplotlib.pyplot as plt
# from tkinter import *



sim = rebound.Simulation.from_file("solar_system.bin")
rebound.OrbitPlot(sim)
# plt.show()

sim.add(m=1,x=40, y=-1000, vy=1.2)
sim.integrate(sim.t+600)

sim.move_to_hel()



for i in range(300):
    sim.integrate(sim.t+10)
    fig = rebound.OrbitPlot(sim,color=True,unitlabel="[AU]", lim=250.,
        show_orbit=True, fancy=True)
    plt.close(fig)
    print(i, "of 300.")
    # plt.savefig(f"plotje_{i}.png")
plt.show()
