


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmasher as cmr

# add directory above to Python Path 
import sys
sys.path.append('../')

import rlightcones as lightcones

model = sys.argv[1]
noise_model = sys.argv[2]


lightcone = lightcones.Lightcone(model, noise_model)

cmap = cmr.torch_r

plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')

left  = 0.15
height = 0.8
bottom = 0.15
width = 0.8

fig = plt.figure(figsize = (3.5, 3.5))

ax = fig.add_axes((left, bottom, width, height))

ax.axhline(5., c='k', ls=':', lw=1, alpha=0.5)

zlim = [0,10]

ax.plot(zlim,zlim,lw=2, c='k',alpha=0.1)
# ax.scatter(lightcone.z, lightcone.pz, s=5)
# ax.hexbin(lightcone.z, lightcone.pz, gridsize=(100,100))
ax.hist2d(lightcone.z, lightcone.pz, bins=(100,100), range = [[0.,10.],[0.,10.]], norm = 'log', cmap = cmap)


ax.set_xlim(zlim)
ax.set_ylim(zlim)

ax.set_xlabel(r'$\rm z $')
ax.set_ylabel(r'$\rm z_{a}$')

fig.savefig(f'figs/z_pz_{model}.pdf')

fig.clf()
