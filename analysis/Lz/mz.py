
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


ref_filter = 'JWST/NIRCam.F200W'
m_range = [23.5, 30.]
dm = 0.25
m_bin_edges = np.arange(*m_range, dm)
m_bin_centres = 0.5*(m_bin_edges[:-1]+m_bin_edges[1:])
redshifts = np.arange(5.5, 10, 1)


fig, axes = plt.subplots(1, 5, figsize=(7, 2), sharey=True, sharex=True)
plt.subplots_adjust(left=0.075, top=0.925, bottom=0.15, right=0.95, wspace=0.0, hspace=0.0)
colours = cmr.take_cmap_colors('cmr.guppy_r', 4, cmap_range=(0.15, 0.85))


z = {}
m = {}

ceers_dir = '/Users/stephenwilkins/Dropbox/Research/data/images/jwst/ceers'
version = '0.2'
subcat = '-high-z.v0.1'
NIRCam_area = 9.2

for pointing in [1, 2, 3, 6]:

    catalogue_id = f'{ceers_dir}/cats/CEERS_NIRCam{pointing}_v{version}{subcat}'
    catalogue_filename = f'{catalogue_id}.h5'
    ref_filter_ = ref_filter.split('.')[-1][:-1]

    with h5py.File(catalogue_filename, 'r') as hf:

        print(pointing, '-'*50)
        hf.visit(print)

        z[pointing] = hf['z'][:]
        m[pointing] = flux_to_m(hf['photom/'+ref_filter_][:])


for ax, redshift in zip(axes.flatten(), redshifts):

    label = rf'$\rm {redshift-0.5:.0f}<z<{redshift+0.5:.0f}$'

    ax.text(0.5, 1.025, label, fontsize=7, color='k', transform=ax.transAxes, ha='center')

    print('-'*20, redshift)

    N = np.zeros(len(m_bin_centres))

    for i, (pointing, ls, colour) in enumerate(zip([1, 2, 3, 6], ['-', '--', '-.', ':'], colours)):

        s = np.fabs(redshift-z[pointing]) < 0.5

        print(pointing, np.sum(s))

        n, _ = np.histogram(m[pointing][s], bins=m_bin_edges)
        N += n
        print(m_bin_centres, N)
        ax.fill_between(m_bin_centres, N*0.0, N, step='mid', color=colour,
                        label=rf'$\rm CEERS/NIRCam{pointing}$', zorder=4-i)

        # ax.step(m_bin_centres, N, where='mid', color=colour, ls=ls,
        #         label=f'NIRCam{pointing}',)

        # ax.plot(m_bin_centres, N, color=colour, ls=ls, lw=1,
        #         label=rf'$\rm NIRCam{pointing}$')

        ax.set_ylim([0, 49.9])

    ax.step(m_bin_centres, N, where='mid',
            ls='-', c='k', lw=1, zorder=5)

axes[0].legend(loc='upper left', fontsize=7, labelspacing=0.1)


axes[0].set_ylabel(r"$\rm N$", fontsize=8)

# ax.set_xlabel(rf"$\rm {ref_filter}$", fontsize=8)

fig.savefig(f'figs/mz.pdf')
