
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cmasher as cmr
import h5py
from astropy.table import Table
from scipy.stats import binned_statistic

plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')


# add directory above to Python Path 
import sys
sys.path.append('../')

import lightcones
import observations


z_range = [4., 10.]
bin_edges = np.arange(*z_range, 0.2)
bin_centres = 0.5*(bin_edges[:-1]+bin_edges[1:])

filters = [f'JWST/NIRCam.{f}' for f in ['F115W','F150W', 'F200W','F277W','F356W','F410M','F444W']]
colours = [[f1, f2] for f1, f2 in zip(filters[:-1], filters[1:])]

# read observations
obs = observations.CEERS()

# intialise figure
fig, axes = plt.subplots(len(filters), 1, figsize = (3.5,6.5), sharex = True)
plt.subplots_adjust(left=0.15, top=0.975, bottom=0.075, right=0.95, wspace=0.0, hspace=0.0)

# define colour colours
colour_colours = cm.rainbow(np.linspace(0, 1, len(filters)-1))

# redshift histogram
ax = axes[0]
ax.hist(obs.z, bins = bin_edges, log = False, alpha = 0.2, color = 'k')
ax.set_ylabel(r'$\rm N$') # colour notation is e.g. Y-J
ax.set_xlim(z_range)


for (f1, f2), ax, color in zip(colours, axes[1:], colour_colours):

    # get array of colours
    c = obs.c(f1, f2)
    
    # calculate statistics

    mmin, _, _ = binned_statistic(obs.z, c, bins = bin_edges, statistic = min)
    mmax, _, _ = binned_statistic(obs.z, c, bins = bin_edges, statistic = max)
    P16, _, _ = binned_statistic(obs.z, c, bins = bin_edges, statistic = lambda x: np.percentile(x, 16))
    P84, _, _  = binned_statistic(obs.z, c, bins = bin_edges, statistic = lambda x: np.percentile(x, 84))
    P5, _, _ = binned_statistic(obs.z, c, bins = bin_edges, statistic = lambda x: np.percentile(x, 5))
    P95, _, _  = binned_statistic(obs.z, c, bins = bin_edges, statistic = lambda x: np.percentile(x, 95))
    med, _, _ = binned_statistic(obs.z, c, bins = bin_edges, statistic = 'median')
    N, _ = np.histogram(obs.z, bins = bin_edges)


    # ax.scatter(obs.z,c,marker='.', s=1, c=color, alpha = 0.5)

    # fill between the P5 and P95 percentiles
    ax.fill_between(bin_centres, P5, P95, color=color, alpha=0.2, lw=0)

    # fill between the P16 and P84 percentiles
    ax.fill_between(bin_centres, P16, P84, color=color, alpha=0.4, lw=0)

    # plot median line
    ax.plot(bin_centres, med, c=color, lw=2, ls='-')

    ax.plot(bin_centres, mmin, c=color, lw=1, ls='-', alpha = 0.5)
    ax.plot(bin_centres, mmax, c=color, lw=1, ls='-', alpha = 0.5)

    ax.set_xlim(z_range)
    ax.set_ylim(-0.99, 1.49)
    ax.set_yticks(np.arange(-0.5, 1.5, 0.5))
    # ax.set_ylabel(rf"$\rm log_{{10}}(f_{{ {f2.split('.')[-1]} }}/f_{{ {f1.split('.')[-1]} }}) $", fontsize=8)


    ax.set_ylabel(rf"$\rm {f1.split('.')[-1]} - {f2.split('.')[-1]}$", fontsize=7)

axes[-1].set_xlabel(rf"$\rm z$", fontsize=8)


fig.savefig(f'figs/cz_obs.pdf')
