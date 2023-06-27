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

plt.rcParams["figure.figsize"] = (3.0, 2.0)

# tableau-colorblind10 (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/tableau-colorblind10.html)
#colorblind_colors = ["#006BA4", "#FF800E", "#ABABAB", "#595959", "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]

# seaborn-colorblind (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/seaborn-colorblind.html)
colorblind_colors = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9"]

particle = pp.particle.GammaDef()
medium = pp.medium.Air()

param_no_lpm = pp.parametrization.photopair.KochMotz()
param_lpm = pp.parametrization.photopair.KochMotz(True, particle, medium)

target = medium.components[0]
assert(target.name == "N")

x_list = np.linspace(0, 1, 10000)

E_list = [1e11, 1e12, 1e13, 1e14, 1e15]
colors = [colorblind_colors[0], colorblind_colors[1], colorblind_colors[2], colorblind_colors[3], colorblind_colors[4], colorblind_colors[5]]

for E, color in zip(E_list, colors):
	logE_pev = np.log10(E/1e9)
	plt.plot(x_list, np.vectorize(param_lpm.differential_crosssection)(particle, target, E, x_list), color=color, label=f'E=\\SI{{e{logE_pev:.0f}}}{{\\peta\\electronvolt}}')

plt.plot(x_list, np.vectorize(param_no_lpm.differential_crosssection)(particle, target, E, x_list), color='tab:red', linestyle='dashed', label='No LPM')


plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2, fontsize=8)
plt.yscale('log')
plt.grid()
plt.ylabel(r'$ \left( \mathrm{d}\sigma \,/\, \mathrm{d} x \right) \,/\, \si{\centi\meter\squared\per\gram} $', fontsize=8)
plt.xlabel(r'$x$', fontsize=8)


plt.savefig('lpm_photopair_differential_small.pdf', bbox_inches='tight')