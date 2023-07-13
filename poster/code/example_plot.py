import proposal as pp

import matplotlib.pyplot as plt
import numpy as np

pp.RandomGenerator.get().set_seed(1234)

# read properties from config file
particle = pp.particle.EMinusDef()
prop = pp.Propagator(particle, "config.json")

# define initial particle state
init_state = pp.particle.ParticleState()
init_state.position = pp.Cartesian3D(0, 0, 0)
init_state.direction = pp.Cartesian3D(0, 0, 1)
init_state.energy = 1e3 # MeV

# propagation
ranges = []
for i in range(int(1e6)):
    output = prop.propagate(init_state)
    x_f = output.final_state().propagated_distance
    ranges.append(x_f)


#bins = np.geomspace(min(ranges), max(ranges), 30)
bins = np.geomspace(1e3, max(ranges), 30)
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.hist(ranges, bins=bins)
plt.title('Electron ranges in air', fontsize=14)
plt.xlabel('propagated distance / cm', fontsize=14)
plt.ylabel('# particles', fontsize=14)
plt.tight_layout()
plt.savefig("example_output.pdf", dpi=300)
