
import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cmasher as cmr
import h5py
from astropy.table import Table
from scipy.stats import binned_statistic

# set style
plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')

# add directory above to Python Path 
import sys
sys.path.append('../')

import rlightcones as lightcones
import observations


z_range = [4., 10.]
bin_edges = np.arange(*z_range, 0.25)
bin_centres = 0.5*(bin_edges[:-1]+bin_edges[1:])




filters = [f'JWST/NIRCam.{f}' for f in ['F115W', 'F150W',
                                        'F200W', 'F277W', 'F356W', 'F410M', 'F444W']]


colours = [[f1, f2] for f1, f2 in zip(filters[:-1], filters[1:])]
fl = {f: f.split('.')[-1] for f in filters}

# get lightcones


models = ['flares', 'scsam', 'jaguar', 'dream']
models = ['scsam','jaguar', 'dream']
noise_model = 'no'
noise_model = 'basic'

models = ['scsam','scsam_10000']


lc = {model: lightcones.Lightcone(model, noise_model) for model in models}

# read observations
obs = observations.CEERS()


fig, axes = plt.subplots(len(filters)-1, 1, figsize=(3.5, 6), sharex=True)
plt.subplots_adjust(left=0.15, top=0.975, bottom=0.05, right=0.95, wspace=0.0, hspace=0.0)


colour_colours = cm.rainbow(np.linspace(0, 1, len(filters)-1))


for i, ((f1, f2), ax, color) in enumerate(zip(colours, axes, colour_colours)):


    # add observations
    
    # get array of colours
    c = obs.c(f1, f2)
    
    # calculate statistics
    med, _, _ = binned_statistic(obs.z, c, bins = bin_edges, statistic = 'median')

    # plot median line
    ax.plot(bin_centres, med, c=color, lw=3, ls='-', alpha = 0.3)

    # add models

    for model, ls in zip(models, ['-', '--', '-.', ':']):

        lightcone = lc[model]

        s = (lightcone.fnu['JWST/NIRCam.F200W'] > 50) 
        print(np.sum(s))

        if f1 in lightcone.m.keys() and f2 in lightcone.m.keys():

            c = lightcone.m[f1] - lightcone.m[f2]

            med, _, _ = binned_statistic(lightcone.z[s], c[s], bins=bin_edges, statistic='median')
            mean, _, _ = binned_statistic(lightcone.z[s], c[s], bins=bin_edges, statistic='mean')
            N, _ = np.histogram(lightcone.z[s], bins=bin_edges)

            ax.plot(bin_centres, med, color=color, lw=1, ls=ls, alpha=0.7, label=rf'$\rm {lightcone.name} $')

            # ax.plot(bin_centres, med, color='k', lw=1, ls=ls,
            #         alpha=0.7, label=rf'$\rm {lightcone.name} (median)$')

            # ax.plot(bin_centres, mean, color='k', lw=2, ls=ls,
            #         alpha=0.7, label=rf'$\rm {lightcone.name} (mean)$')

    if i == 0:
        ax.legend(fontsize=7, labelspacing=0.1)

    ax.set_xlim(z_range)
    ax.set_ylim(-0.99, 1.49)
    ax.set_yticks(np.arange(-0.5, 1.5, 0.5))
    # ax.set_ylabel(rf"$\rm log_{{10}}(f_{{ {f2.split('.')[-1]} }}/f_{{ {f1.split('.')[-1]} }}) $", fontsize=8)

    ax.set_ylabel(rf"$\rm {fl[f1]} - {fl[f2]}$", fontsize=7)

axes[-1].set_xlabel(rf"$\rm z$", fontsize=8)

fig.savefig(f'figs/cz_models_{noise_model}noise.pdf')
