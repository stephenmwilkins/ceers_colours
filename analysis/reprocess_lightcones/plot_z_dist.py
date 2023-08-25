


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmasher as cmr

# add directory above to Python Path 
import sys
sys.path.append('../')

import rlightcones as lightcones

model = 'scsam'
noise_model = 'basic'
lightcone = lightcones.Lightcone(model, noise_model)

s = lightcone.get_selection()

cmap = cmr.torch_r

plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')

left  = 0.15
height = 0.8
bottom = 0.15
width = 0.8

fig = plt.figure(figsize = (3.5, 2.5))

ax = fig.add_axes((left, bottom, width, height))

zlim = [0,10]


# true
ax.hist(lightcone.z[s], bins = 100, range = [0,10.], log=True, histtype = 'stepfilled', color='k', alpha = 0.3, label = r'$\rm true\ z$')

# photo-z
ax.hist(lightcone.pz[s], bins = 100, range = [0,10.], log=True, histtype = 'step', color='k', label = r'$\rm photo-z$')

ax.legend(fontsize=8)
ax.set_xlim(zlim)

ax.set_xlabel(r'$\rm z $')
ax.set_ylabel(r'$\rm N$')

fig.savefig(f'figs/z_dist_{model}.pdf')

fig.clf()
