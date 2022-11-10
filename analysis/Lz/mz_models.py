
import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import cmasher as cmr

import h5py
from astropy.table import Table

from flare.photom import flux_to_m
import flares
import flares_analysis
import flare.plt as fplt
from scipy.stats import binned_statistic

from flare_lf import lightcones


ref_filter = 'JWST/NIRCam.F200W'
m_range = [23.5, 28.5]
dm = 0.25
m_bin_edges = np.arange(*m_range, dm)
m_bin_centres = 0.5*(m_bin_edges[:-1]+m_bin_edges[1:])
redshifts = np.arange(5.5, 10, 1)


models = ['scsam', 'jaguar', 'dream']  # 'flares',
# models = ['scsam']


lc = {model: lightcones.Lightcone(model) for model in models}


fig, axes = plt.subplots(1, 5, figsize=(7, 2), sharey=True, sharex=True)
plt.subplots_adjust(left=0.075, top=0.925, bottom=0.15, right=0.95, wspace=0.0, hspace=0.0)
colours = cmr.take_cmap_colors('cmr.guppy_r', len(redshifts), cmap_range=(0.15, 0.85))


# ceers = ReadCEERS()
# s = ceers.get_selection()
# z = ceers.z[s]
# m = ceers.m[ref_filter.split('.')[-1]][s]

ceers_dir = '/Users/stephenwilkins/Dropbox/Research/data/images/jwst/ceers'

version = '0.2'
subcat = '-high-z.v0.1'

z = np.array([])
m = np.array([])

for pointing in [1, 2, 3, 6]:

    catalogue_id = f'{ceers_dir}/cats/CEERS_NIRCam{pointing}_v{version}{subcat}'
    catalogue_filename = f'{catalogue_id}.h5'
    ref_filter_ = ref_filter.split('.')[-1][:-1]

    with h5py.File(catalogue_filename, 'r') as hf:
        z = np.hstack((z, hf['z'][:]))
        m = np.hstack((m, flux_to_m(hf['photom/'+ref_filter_][:])))


for ax, redshift, colour in zip(axes.flatten(), redshifts, colours):

    label = rf'$\rm {redshift-0.5:.0f}<z<{redshift+0.5:.0f}$'

    ax.text(0.5, 1.025, label, fontsize=7, color=colour, transform=ax.transAxes, ha='center')

    # --- field by field observations
    s_ = np.fabs(redshift-z) < 0.5
    N, _ = np.histogram(m[s_], bins=m_bin_edges)
    ax.step(m_bin_centres, N/32., where='mid', color=colour)

    # --- summed observtions

    # --- models

    for model, ls in zip(models, ['-', '--', '-.', ':']):

        lightcone = lc[model]

        s_ = np.fabs(redshift-lightcone.z) < 0.5

        N, _ = np.histogram(lightcone.m[ref_filter][s_], bins=m_bin_edges)

        ax.plot(m_bin_centres, N/lightcone.area, color='k', lw=1,
                ls=ls, alpha=0.7, label=rf'$\rm {lightcone.name}$')

    # ax.set_xlim(z_range)
    # ax.set_ylim(-0.99, 1.49)

    ax.set_yscale('log')

axes[-1].legend(fontsize=7, labelspacing=0.1)
axes[0].set_ylabel(r"$\rm N/arcmin^{2}$", fontsize=7)
# ax.set_xlabel(rf"$\rm {ref_filter}$", fontsize=8)
#

# ax.legend(fontsize=7, labelspacing = 0.1)
# axes[-1].set_xlabel(rf"$\rm z$", fontsize=8)

fig.savefig(f'figs/mz_models.pdf')
