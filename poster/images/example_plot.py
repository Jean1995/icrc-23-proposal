import proposal as pp

import matplotlib.pyplot as plt
import numpy as np

pp.RandomGenerator.get().set_seed(1234)

# read properties from config file
particle = pp.particle.MuMinusDef()
prop = pp.Propagator(particle, "config.json")

# define initial particle state
init_state = pp.particle.ParticleState()
init_state.position = pp.Cartesian3D(0, 0, 0)
init_state.direction = pp.Cartesian3D(0, 0, 1)
init_state.energy = 1e9 # MeV

# propagation
final_energies = []
for i in range(10000):
    output = prop.propagate(init_state, 
                            max_distance = 1e5) # cm
    E_f = output.final_state().energy
    final_energies.append(E_f)


bins = np.geomspace(min(final_energies), max(final_energies), 50)
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.hist(final_energies, bins=bins)
plt.title('Muon energies after 1 km of ice')
plt.xlabel('final_energies / MeV')
plt.ylabel('# particles')
plt.tight_layout()
plt.savefig("example_output.pdf", dpi=300)
