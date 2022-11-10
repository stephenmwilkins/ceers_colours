
import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import cmasher as cmr

import h5py
from astropy.table import Table

import flares
import flares_analysis
import flare.plt as fplt
from scipy.stats import binned_statistic

from flare_lf import lightcones


z_range = [5., 10.]
bin_edges = np.arange(*z_range, 0.1)
bin_centres = 0.5*(bin_edges[:-1]+bin_edges[1:])

filters = [f'JWST/NIRCam.{filter}' for filter in ['F115W','F150W', 'F200W','F277W','F356W','F410M','F444W']]
colours = [[f1, f2] for f1, f2 in zip(filters[:-1], filters[1:])]
fl = {f: f.split('.')[-1] for f in filters}

# --- open observations

template_set = 'Larson22'


models = ['jaguar']

for model in models:


    hf = h5py.File(f'data/{model}.hf', 'r')

    z = hf['input/z'][()]
    pz = hf[f'pz/eazy/{template_set}/z_a'][()]

    fig, axes = plt.subplots(len(filters), 1, figsize = (3.5,6), sharex = True)
    plt.subplots_adjust(left=0.15, top=0.975, bottom=0.05, right=0.95, wspace=0.0, hspace=0.0)


    colour_colours = cm.rainbow(np.linspace(0, 1, len(filters)-1))

    # --- redshift histogram
    ax = axes[0]
    ax.hist(z, bins = bin_edges, log = False, alpha = 0.2, color = 'k')
    ax.set_ylabel(r'$\rm N$') # colour notation is e.g. Y-J
    ax.set_xlim(z_range)

    for i, ((f1, f2), ax, color) in enumerate(zip(colours, axes[1:], colour_colours)):

        c = -2.5*np.log10(hf[f'input/{f2}/flux'][()]/hf[f'input/{f1}/flux'][()])

        med, _, _ = binned_statistic(z, c, bins = bin_edges, statistic = 'median')
        N, _ = np.histogram(z, bins = bin_edges)

        ax.plot(bin_centres, med, color='k', lw=1, ls='-', alpha = 0.7)

        c = -2.5*np.log10(hf[f'obs/{f2}/flux'][()]/hf[f'obs/{f1}/flux'][()])

        med, _, _ = binned_statistic(pz, c, bins = bin_edges, statistic = 'median')
        N, _ = np.histogram(z, bins = bin_edges)

        ax.plot(bin_centres, med, color='k', lw=1, ls='--', alpha = 0.7)


        ax.set_xlim(z_range)
        ax.set_ylim(-0.99, 1.49)
        ax.set_yticks(np.arange(-0.5, 1.5, 0.5))
        # ax.set_ylabel(rf"$\rm log_{{10}}(f_{{ {f2.split('.')[-1]} }}/f_{{ {f1.split('.')[-1]} }}) $", fontsize=8)

        ax.set_ylabel(rf"$\rm {fl[f1]} - {fl[f2]}$", fontsize=8)


    fig.savefig(f'figs/cz_{model}.pdf')
