


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmasher as cmr
from scipy import stats

# add directory above to Python Path 
import sys
sys.path.append('../')

import rlightcones as lightcones

model = sys.argv[1]

lightcones = [lightcones.Lightcone(model, noise_model) for noise_model in ['basic', 'linear']]

filter = 'JWST/NIRCam.F277W'

plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')

left  = 0.15
height = 0.8
bottom = 0.15
width = 0.8

fig = plt.figure(figsize = (3.5, 3.5))

ax = fig.add_axes((left, bottom, width, height))

flux_lims = np.array([0,8])

# zlim = [0,10]

# ax.plot(zlim,zlim,lw=2, c='k',alpha=0.1)
# ax.scatter(np.log10(lightcone.fnu[filter]), np.log10(lightcone.fnu[filter]/lightcone.fnu_err[filter]), s=1)


for lightcone, cmap in zip(lightcones, [cmr.torch_r, cmr.horizon_r]):

    ax.hist2d(np.log10(lightcone.fnu[filter]), np.log10(lightcone.fnu[filter]/lightcone.fnu_err[filter]), bins=(100,100), range = [[0.,10.],[0.,10.]], norm = 'log', cmap = cmap, alpha = 0.5)

    res = stats.linregress(np.log10(lightcone.fnu[filter]), np.log10(lightcone.fnu[filter]/lightcone.fnu_err[filter]))

    ax.plot(flux_lims, res.intercept + res.slope*flux_lims, label=lightcone.noise_model, c=cmap(0.5))



ax.legend(fontsize = 8)

ax.set_xlim([1.,8.])
ax.set_ylim([1.,8.])

ax.set_xlabel(rf'$\rm log_{{10}}(flux({filter})/nJy)$')
ax.set_ylabel(r'$\rm log_{10}(SNR)$')

fig.savefig(f'figs/snr_{model}.pdf')
fig.clf()
