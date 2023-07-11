import proposal as pp
from packaging import version


### check correct proposal version. uncomment this line if you want to run the scripts with newer/other versions
assert version.parse(pp.__version__).major == 7, "You need to install PROPOSAL 7 to run this script"
assert version.parse(pp.__version__).minor >= 6, "You need to install PROPOSAL 7.6 or newer to run this script"
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc

import scipy.constants as c

### conversion factor from mass attenuation coefficient to cross section
# mu/rho = sigma_tot / (u * A)
u = c.u * 1000 # unit: g
A = 28.97 # g/mol
cm2_to_barn = 1e24
conversion = u * A * cm2_to_barn

rc('font', **{'family': 'serif',
   'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('text.latex', preamble="\n".join([r'\usepackage{calrsfs}', r"\usepackage[detect-all,locale=US]{siunitx}"]))

import pandas as pd

# tableau-colorblind10 (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/tableau-colorblind10.html)
#colorblind_colors = ["#006BA4", "#FF800E", "#ABABAB", "#595959", "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]

# seaborn-colorblind (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/seaborn-colorblind.html)
colorblind_colors = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9"]

# plots-specific size (width should always be taken from conf)
plt.rcParams["figure.figsize"] = (6, 3.3)


particle = pp.particle.GammaDef()
medium = pp.medium.Air()
cuts = None

air_nist = pd.read_table("nist_tables/air.txt", sep='  ', header=None, names=["energy", "mu", "mu_en"])

stdcross = pp.crosssection.make_std_crosssection(particle, medium, None, False)


fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.1})

labels = ['Photoeffect', 'Compton', 'Photopair']
cross = [stdcross[3], stdcross[1], stdcross[0] ]
colors = [colorblind_colors[0], colorblind_colors[1], colorblind_colors[2]]

def total_dNdx(E):
	dNdx = 0
	for c in cross:
		dNdx += c.calculate_dNdx(E)
	return dNdx

energies = np.geomspace(min(air_nist['energy']), max(air_nist['energy']), 1000)

for c, label, color in zip(cross, labels, colors):
	ax1.plot(energies, c.calculate_dNdx(energies) * conversion, label=label, color=color)


ax1.plot(energies, total_dNdx(energies) * conversion, label=r'$\sigma_\text{tot}$', color=colorblind_colors[3])


ax1.plot(air_nist['energy'], air_nist["mu"] * conversion, 'x', label='NIST data', color='k', markersize=4)
ax2.plot(air_nist['energy'], air_nist["mu"] / total_dNdx(air_nist['energy']), 'x', label=r'NIST data / $\sigma_\text{tot}$', color='k', markersize=4)


#ax1.set_xlim(3e-3, 4e-3)

#ax1.legend(loc='best')

handles,labels = ax1.get_legend_handles_labels()

#handles = [handles[0], handles[4], handles[1], handles[3], handles[2]]
#labels = [labels[0], labels[4], labels[1], labels[3], labels[2]]

ax1.legend(handles, labels, bbox_to_anchor=(0, 1.01, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, fontsize=10)
ax2.legend(loc='best')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid()
ax1.set_ylabel(r'$\sigma \,/\, \si{\barn} $', fontsize=10)

ax2.grid()
ax2.set_ylim(0.8, 1.2)
ax2.set_ylabel('ratio')

ax2.set_xlabel(r'$E \,/\, \si{\mega\electronvolt}$', fontsize=10)


plt.savefig('photoeffect_nist.pdf', bbox_inches='tight')