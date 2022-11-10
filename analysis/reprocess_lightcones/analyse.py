
from pathlib import Path

import numpy as np
import h5py
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import cmasher as cmr
from scipy.stats import binned_statistic, binned_statistic_2d

from flags.pz import eazy
import flare.plt as fplt
# from flare.filters import add_filters




# --- open h5py catalogue


def simple_plt(figsize = (3.5, 3.5), left = 0.15, bottom = 0.15, width = 0.8, height = 0.8):

    fig = plt.figure(figsize = figsize)
    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax






class Analyser:


    def __init__(self, model, template_set = 'tweak_fsps_QSF_12_v3'):

        print('-'*120)

        self.model = model

        self.cat = h5py.File(f'data/{model}.hf', 'r')


        self.z = self.cat['input/z'][()]
        self.pz = self.cat[f'pz/eazy/{template_set}/z_a'][()]
        self.dz = (self.z - self.pz)/(1+self.z)


    def explore(self):
        self.cat.visit(print) # explore datasets

    def print_info(self, i):

        def visitor_func(name, node):
            if isinstance(node, h5py.Dataset):
                 print(name, node[i])

        self.icat.visititems(visitor_func)


    def dz_plot(self, tag = ''):

        fig, ax = simple_plt(figsize = (3.5, 2.5))


        range_x = [5, 10]
        range_y = [-0.99, 0.99]


        hist, bin_edges_x, bin_edges_y = np.histogram2d(self.z, self.dz, range = [range_x,range_y], bins = (100,100))

        ax.imshow(hist.T, origin = 'lower', aspect = 'auto', extent = [*range_x, *range_y], cmap = cmr.sunburst_r)

        # ax.hist2d(self.z, self.dz, range = [[7, 17],[-1.5, 1.5]], bins = (100,100))

        ax.axhline(0.0, c='k', lw=2, alpha = 0.05)

        for x in [0.1, -0.1]:
            ax.axhline(x, c='k', lw=1, alpha = 0.2, ls = ':')

        ax.set_xlabel(r'$\rm z $')
        ax.set_ylabel(r'$\rm dz/(1+z)$')

        ax.set_yticks([-0.5, 0.0, 0.5])

        fn = f'figs/dz_{self.model}.pdf'
        print(fn)
        fig.savefig(fn)












if __name__ == "__main__":


    model = f'jaguar'

    a = Analyser(model, template_set = 'Larson22')
    a.dz_plot()
