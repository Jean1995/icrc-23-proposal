import proposal as pp
from packaging import version


### check correct proposal version. uncomment this line if you want to run the scripts with newer/other versions
assert version.parse(pp.__version__).major == 7, "You need to install PROPOSAL 7 to run this script"
assert version.parse(pp.__version__).minor >= 6, "You need to install PROPOSAL 7.6 or newer to run this script"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib import rc
rc('font', **{'family': 'serif',
   'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('text.latex', preamble="\n".join([r'\usepackage{calrsfs}', r"\usepackage[detect-all,locale=US]{siunitx}"]))

plt.rcParams["figure.figsize"] = (6, 3.0)

# tableau-colorblind10 (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/tableau-colorblind10.html)
#colorblind_colors = ["#006BA4", "#FF800E", "#ABABAB", "#595959", "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]

# seaborn-colorblind (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/seaborn-colorblind.html)
colorblind_colors = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9"]


photon = pp.particle.GammaDef()
medium = pp.medium.Air()
cuts = None

# from https://www.pdas.com/atmosTable2SI.html
rho_list = np.array([0.089, 0.414, 0.736, 1.225]) / 1000 # in g / cm^3
height_list = np.array([20, 10, 5, 0]) # in km

energies = np.geomspace(1e10, 1e14)

param_photopair_no_lpm = pp.parametrization.photopair.KochMotz(False)
cross_photopair_no_lpm = pp.crosssection.make_crosssection(param_photopair_no_lpm, photon, medium, cuts, False)

conversion_energy = 1e9 # plot energies in PeV

for rho, height, color in zip(rho_list, height_list, colorblind_colors):
    density_correction = rho / rho_list[0] # calculate density correction
    param_photopair_lpm = pp.parametrization.photopair.KochMotz(True, photon, medium, density_correction)
    cross_photopair_lpm = pp.crosssection.make_crosssection(param_photopair_lpm, photon, medium, cuts, False)
    plt.plot(energies/conversion_energy, energies * cross_photopair_lpm.calculate_dNdx(energies), color=color, label=f'$h=\\SI{{{height}}}{{\\kilo\\meter}} \, ( \\rho = \\SI{{{1000*rho}}}{{\\kilo\\gram\\per\\meter\\cubed}})$')

plt.plot(energies/conversion_energy, energies * cross_photopair_no_lpm.calculate_dNdx(energies), color=colorblind_colors[-1], label=r'No LPM')

plt.legend(loc='best')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.ylabel(r'$\sigma \cdot E \,/\, \si{\mega\electronvolt\centi\meter\squared\per\gram} $')
plt.xlabel(r'$E / \si{\peta\electronvolt}$')


plt.savefig('lpm_cross_photopair.pdf', bbox_inches='tight')