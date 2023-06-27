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

# plots-specific size (width should always be taken from conf)
plt.rcParams["figure.figsize"] = (6, 2.5)


# tableau-colorblind10 (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/tableau-colorblind10.html)
#colorblind_colors = ["#006BA4", "#FF800E", "#ABABAB", "#595959", "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]

# seaborn-colorblind (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/seaborn-colorblind.html)
colorblind_colors = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9"]

particle = pp.particle.GammaDef()
medium = pp.medium.Air()
cuts = None

param_list = []

param_list.append(pp.parametrization.photoproduction.BezrukovBugaev())
param_list.append(pp.parametrization.photoproduction.Caldwell())
param_list.append(pp.parametrization.photoproduction.Kokoulin())
param_list.append(pp.parametrization.photoproduction.Zeus())
param_list.append(pp.parametrization.photoproduction.Rhode())
param_list.append(pp.parametrization.photoproduction.Heck())
param_list.append(pp.parametrization.photoproduction.HeckC7Shadowing())

labels = [r"Bezrukov \& Bugaev", "Caldwell", "Kokoulin", "Zeus", "Rhode", "Heck"]
colors = colorblind_colors

cross_list = []
for param in param_list:
	cross_list.append(pp.crosssection.make_crosssection(param, particle, medium, None, False))

energies = np.geomspace(1.5e2, 1e11, 500)

for cross, label, color in zip(cross_list[:-1], labels, colors):
	plt.plot(energies, cross.calculate_dNdx(energies), label=label, color=color, lw=1)

# with CORSIKA 7 shadowing
plt.plot(energies, cross_list[-1].calculate_dNdx(energies), color=colors[-1], lw=1, linestyle='dashed')

plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3)
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.ylabel(r'$\sigma \,/\, \si{\centi\meter\squared\per\gram} $')
plt.xlabel(r'$E \,/\, \si{\mega\electronvolt}$')
plt.ylim(3e-5, 3e-4)
#plt.xlim(1.5e2, 1e11)


plt.savefig('photoproduction_cross.pdf', bbox_inches='tight')