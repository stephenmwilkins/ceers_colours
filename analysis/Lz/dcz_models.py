
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

models = ['flares', 'scsam', 'jaguar', 'dream']
models = ['flares']

filters = [f'JWST/NIRCam.{f}' for f in ['F115W','F150W', 'F200W','F277W','F356W','F410M','F444W']]
colours = [[f1, f2] for f1, f2 in zip(filters[:-1], filters[1:])]
fl = {f: f.split('.')[-1] for f in filters}

# --- open observations

from ceers import ReadCEERS
d = ReadCEERS()
s = d.get_selection()
z = d.z[s]


lc = {model: lightcones.Lightcone(model) for model in ['flares', 'scsam', 'jaguar', 'dream']}




fig, axes = plt.subplots(len(filters)-1, 1, figsize = (3.5,6), sharex = True)
plt.subplots_adjust(left=0.15, top=0.975, bottom=0.05, right=0.95, wspace=0.0, hspace=0.0)

colour_colours = cm.rainbow(np.linspace(0, 1, len(filters)-1))


for i, ((f1, f2), ax, color) in enumerate(zip(colours, axes, colour_colours)):

    # --- observations
    f1_ = f1.split('.')[-1]
    f2_ = f2.split('.')[-1]

    c = d.m[f1_][s] - d.m[f2_][s]


    obs, _, _ = binned_statistic(z, c, bins = bin_edges, statistic = 'median')


    for model, ls in zip(models, ['-','--','-.',':']):

        lightcone = lc[model]

        if f1 in lightcone.m.keys() and f2 in lightcone.m.keys():

            c = lightcone.m[f1] - lightcone.m[f2]

            mod, _, _ = binned_statistic(lightcone.z, c, bins = bin_edges, statistic = 'median')

            print(mod-obs)

            ax.plot(bin_centres, mod - obs, color=color, lw=1, ls=ls, alpha = 0.7, label = rf'$\rm {lightcone.name}$')


    if i==0:
        ax.legend(fontsize=7)

    ax.axhline(0.0, lw=3, c='k', alpha = 0.1)
    ax.set_xlim(z_range)
    ax.set_ylim(-0.49, 0.49)
    # ax.set_yticks(np.arange(-0.5, 1.5, 0.5))
    # ax.set_ylabel(rf"$\rm log_{{10}}(f_{{ {f2.split('.')[-1]} }}/f_{{ {f1.split('.')[-1]} }}) $", fontsize=8)


    ax.set_ylabel(rf"$\rm {fl[f1]} - {fl[f2]}$", fontsize=8)


fig.savefig(f'figs/dcz_models.pdf')
