
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

import observations



ref_filter = 'JWST/NIRCam.F277W'
m_range = [23.5, 30.]
dm = 0.25
m_bin_edges = np.arange(*m_range, dm)
m_bin_centres = 0.5*(m_bin_edges[:-1]+m_bin_edges[1:])
redshifts = np.arange(5.5, 10, 1)


fig, axes = plt.subplots(1, 5, figsize=(7, 2), sharey=True, sharex=True)
plt.subplots_adjust(left=0.075, top=0.925, bottom=0.15, right=0.95, wspace=0.0, hspace=0.0)
colours = cmr.take_cmap_colors('cmr.guppy_r', 4, cmap_range=(0.15, 0.85))


# read observations
obs = observations.CEERS()




for ax, redshift in zip(axes.flatten(), redshifts):

    label = rf'$\rm {redshift-0.5:.0f}<z<{redshift+0.5:.0f}$'

    ax.text(0.5, 1.025, label, fontsize=7, color='k', transform=ax.transAxes, ha='center')

    print('-'*20, redshift)

    # select galaxies from the observations
    s = np.fabs(obs.z-redshift)<0.5

    m = obs.m[ref_filter]

    N, _ = np.histogram(m[s], bins=m_bin_edges)
  
    ax.set_ylim([0, 49.9])

    ax.step(m_bin_centres, N, where='mid',
            ls='-', c='k', lw=1, zorder=5)

axes[0].legend(loc='upper left', fontsize=7, labelspacing=0.1)

axes[0].set_ylabel(r"$\rm N$", fontsize=8)

# ax.set_xlabel(rf"$\rm {ref_filter}$", fontsize=8)

fig.savefig(f'figs/mz_obs.pdf')
