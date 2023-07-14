import matplotlib.pyplot as plt
import numpy as np
import proposal as pp

import pandas as pd
import scipy.constants as c

from matplotlib import rc
rc('font', **{'family': 'serif',
   'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('text.latex', preamble="\n".join([r'\usepackage{calrsfs}', r"\usepackage[detect-all,locale=US]{siunitx}"]))

plt.rcParams["figure.figsize"] = (6, 2.9)

# tableau-colorblind10 (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/tableau-colorblind10.html)
#colorblind_colors = ["#006BA4", "#FF800E", "#ABABAB", "#595959", "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]

# seaborn-colorblind (https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/seaborn-colorblind.html)
colorblind_colors = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9"]


def calculate_losses_photon(energies, medium, ecut=np.inf, modus='dedx'):

    p_def = pp.particle.GammaDef()
    if ecut == 0:
        ecuts = None
    else:
        ecuts = pp.EnergyCutSettings(ecut, 1, False)

    args = {
        "particle_def": p_def,
        "target": medium,
        "interpolate": False,
        "cuts": None,
    }
    lpm_args = {
        "particle_def": p_def,
        "lpm": True,
        "medium": medium,
    }
    nolpm_args = {
        "particle_def": p_def,
        "lpm": False,
        "medium": medium,
    }
    args_ecut =  {**args, **{"cuts": ecuts}}

    #   Create cross sections out of their parametrizations
    cross_photopair = pp.crosssection.make_crosssection(
        parametrization=pp.parametrization.photopair.KochMotz(**lpm_args),
        **args)
    cross_photopair_nolpm = pp.crosssection.make_crosssection(
        parametrization=pp.parametrization.photopair.KochMotz(**nolpm_args),
        **args)
    cross_compton = pp.crosssection.make_crosssection(
        parametrization=pp.parametrization.compton.KleinNishina(),
        **args_ecut)
    cross_photoprod = pp.crosssection.make_crosssection(
        parametrization=pp.parametrization.photoproduction.Rhode(),
        **args)
    cross_photoeffect = pp.crosssection.make_crosssection(
        parametrization=pp.parametrization.photoeffect.Sauter(),
        **args)
    cross_photomupair = pp.crosssection.make_crosssection(
        parametrization=pp.parametrization.photomupair.BurkhardtKelnerKokoulin(),
        **args)
    cross_list = [cross_photoeffect, cross_compton, cross_photopair, cross_photopair_nolpm, cross_photoprod, cross_photomupair]

    #   Calculate dE/dx or dN/dx at the given energies
    xsection = np.empty((len(cross_list), len(energies)))

    for idx, cross in enumerate(cross_list):
        if modus == 'dedx':
            if cross.param_name == 'KleinNishina':
                xsection[idx] = cross.calculate_dEdx(energies)
            else:
                xsection[idx] = energies * cross.calculate_dNdx(energies)
        elif modus == 'dndx':
            xsection[idx] = cross.calculate_dNdx(energies)
        else:
            raise AttributeError('modus must be dedx or dndx')

    return xsection




def plot_dEdx(energy_min=1e3, energy_max=1e12, n_energies=200, ecut=np.inf,
    modus='dedx', bottom=None, particle='Muon', medium='Ice'):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    energies = np.logspace(
        np.log10(energy_min),
        np.log10(energy_max),
        n_energies)

    ### conversion factor from mass attenuation coefficient to cross section
    # mu/rho = sigma_tot / (u * A)
    u = c.u * 1000 # unit: g
    cm2_to_barn = 1e24

    if medium == 'Ice':
        medium_proposal = pp.medium.Ice()
        A = 267.02047 # g/mol
    elif medium == 'Air':
        medium_proposal = pp.medium.Air()
        A = 28.97 # g/mol
    else:
        raise AttributeError('medium must be medium or air (or you need to add it to the script)')

    if modus == 'dndx':
        conversion = u * A * cm2_to_barn
    else:
        conversion = 1


    loss_sum = np.zeros(len(energies))

    if particle == 'Photon':
        dedx_all = calculate_losses_photon(energies, medium_proposal, ecut, modus=modus)
        labels = [
            'Photoelectric absorption',
            'Compton',
            r'$e^+ e^-$ pair production (LPM on)',
            r'$e^+ e^-$ pair production (LPM off)',
            'Photonuclear interactions',
            r'$\mu^+ \mu^-$ pair production',
        ]
        colors = [colorblind_colors[0], colorblind_colors[1], colorblind_colors[2], colorblind_colors[2],
         colorblind_colors[4], colorblind_colors[3]]

    else:
        print("Not implemented...")

    for dedx, _label, _color in zip(dedx_all, labels, colors):
        if "LPM off" in _label:
            _linestyle = '--'
        else:
            _linestyle = '-'
            loss_sum += dedx
        ax.plot(energies, dedx * conversion, linestyle=_linestyle, label=_label, color=_color)

    #if (modus == 'dndx' and particle=='Photon'):
    #    #plot nist tables
    #    air_nist = pd.read_table("nist_tables/air.txt", sep='  ', header=None, names=["energy", "mu", "mu_en"])
    #    ax.plot(air_nist['energy'], air_nist["mu"],linestyle='dashed', label='NIST data', color=colorblind_colors[5])


    # use this if you want to plot the sum of all cross sections
    # ax.plot(energies, loss_sum, linestyle='-', label="Sum", color="tab:gray")

    if modus == 'dedx':
        output ='{particle}_{medium}_{modus}.pdf'
        ax.set_ylabel(r'$\left\langle -\frac{\mathrm{d}E}{\mathrm{d}X}\right\rangle \,\left/\, \left( \mathrm{MeV} \mathrm{g}^{-1} \mathrm{cm}^2 \right) \right. $')
    elif modus == 'dndx':
        output ='{particle}_{medium}_{modus}_ecut_{ecut:.4g}.pdf'
        ax.set_ylabel(r'$\sigma  \,/\, \mathrm{b} $', fontsize=10)

    output = output.format(**{'ecut':ecut, 'modus':modus, 'particle':particle, 'medium':medium})
    if bottom:
        ax.set_ylim(bottom=bottom)
    ax.set_xlabel(r'$E \,/\, \mathrm{MeV} $', fontsize=10)
    #ax.legend(loc='best')

    handles,labels = ax.get_legend_handles_labels()

    handles = [handles[0], handles[3], handles[2], handles[4], handles[1], handles[5]]
    labels = [labels[0], labels[3], labels[2], labels[4], labels[1], labels[5]]


    ax.legend(handles, labels, bbox_to_anchor=(0, 1.01, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2, fontsize=10)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(min(energies), max(energies))
    plt.grid()
    plt.minorticks_on()
    plt.savefig(output, bbox_inches='tight', pad_inches=0.02)
    plt.close()



if __name__ == "__main__":
    plot_dEdx(ecut=0, modus='dndx', particle='Photon', bottom=1e-8, energy_min=1e-3, medium='Air', energy_max=1e14)