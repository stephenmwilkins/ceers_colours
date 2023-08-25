
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


ref_filter = 'JWST/NIRCam.F277W'
m_range = [23.5, 28.5]
dm = 0.25
m_bin_edges = np.arange(*m_range, dm)
m_bin_centres = 0.5*(m_bin_edges[:-1]+m_bin_edges[1:])
redshifts = np.arange(5.5, 10, 1)

noise_model = 'no'
# noise_model = 'basic'
models = ['scsam', 'jaguar', 'dream']  # 'flares',
# models = ['scsam']

# read observations
obs = observations.CEERS()




# get lightcones
lc = {model: lightcones.Lightcone(model, noise_model) for model in models}

fig, axes = plt.subplots(1, 5, figsize=(7, 2), sharey=True, sharex=True)
plt.subplots_adjust(left=0.075, top=0.925, bottom=0.2, right=0.95, wspace=0.0, hspace=0.0)
colours = cmr.take_cmap_colors('cmr.guppy_r', len(redshifts), cmap_range=(0.15, 0.85))



for ax, redshift, colour in zip(axes.flatten(), redshifts, colours):

    label = rf'$\rm {redshift-0.5:.0f}<z<{redshift+0.5:.0f}$'

    ax.text(0.5, 1.025, label, fontsize=7, color=colour, transform=ax.transAxes, ha='center')

    # add models

    for model, ls in zip(models, ['--', '-.', ':']):

        lightcone = lc[model]

        # select galaxies the same as the observations
        s = lightcone.get_selection()

        s = s & (np.fabs(redshift-lightcone.z) < 0.5)
        
        N, _ = np.histogram(lightcone.m[ref_filter][s], bins=m_bin_edges)

        ax.plot(m_bin_centres, N/lightcone.area/dm, color=colour, lw=1,
                ls=ls, alpha=0.7, label=rf'$\rm {lightcone.name}$')


    # add observations

    s = np.fabs(obs.z-redshift)<0.5

    m = obs.m[ref_filter]

    N, _ = np.histogram(m[s], bins=m_bin_edges)
  
    ax.step(m_bin_centres, N/obs.area/dm, where='mid',
            ls='-', c=colour, lw=2, alpha=0.5, zorder=5)
  
    ax.axhline(1/obs.area/dm,c='k',lw=1, alpha=0.2)
    ax.set_ylim(0.01,20)
    ax.set_yscale('log')

    

axes[-1].legend(fontsize=7, labelspacing=0.1)
axes[0].set_ylabel(r"$\rm N/arcmin^{2}\ mag^{-1}$", fontsize=7)
axes[2].set_xlabel(rf"$\rm F277W$", fontsize=8)

fig.savefig(f'figs/mz_models_{noise_model}noise.pdf')
