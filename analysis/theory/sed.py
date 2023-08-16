
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmasher as cmr

from synthesizer.filters import FilterCollection
from synthesizer.grid import Grid
from synthesizer.parametric.sfzh import SFH, ZH, generate_sfzh
from synthesizer.parametric.galaxy import Galaxy
from synthesizer.plt import single, single_histxy, mlabel
from unyt import yr, Myr
from synthesizer.igm import Madau96, Inoue14
from astropy.cosmology import Planck18 as cosmo

plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')



# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)

grid_name = "bpass-2.2.1-bin_chabrier03-0.1,300.0_cloudy-c17.03"
grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids'
grid = Grid(grid_name, grid_dir=grid_dir)


# -------------------------------------------------
# --- define star formation and metal enrichment history (sfzh)

# define filters
filter_codes = [f'JWST/NIRCam.{f}' for f in ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']]  # define a list of filter codes
fc = FilterCollection(filter_codes, new_lam=grid.lam)

print(fc.filters)

# define the parameters of the star formation and metal enrichment histories
sfh_p = {'duration': 100 * Myr}
Z_p = {'log10Z': -2.0}  # can also use linear metallicity e.g. {'Z': 0.01}
stellar_mass = 1E8

# define the functional form of the star formation and metal enrichment histories
sfh = SFH.Constant(sfh_p)  # constant star formation
Zh = ZH.deltaConstant(Z_p)  # constant metallicity

# get the 2D star formation and metal enrichment history for the given SPS grid. This is (age, Z).
sfzh = generate_sfzh(grid.log10age, grid.metallicity, sfh, Zh, stellar_mass=stellar_mass)

# create a galaxy object
galaxy = Galaxy(sfzh)

# generate spectra using pacman model (complex)
sed = galaxy.get_spectra_intrinsic(grid)




colours = cm.rainbow(np.linspace(0, 1, len(fc.filters)))
colour_colours = cm.rainbow(np.linspace(0, 1, len(fc.filters)-1))


left  = 0.15
height = 0.75
bottom = 0.075
width = 0.8
hheight = 0.15

fig = plt.figure(figsize = (3.5, 5.5))

ax = fig.add_axes((left, hheight+bottom, width, height))
filt_ax = fig.add_axes([left, bottom, width, hheight])


for z, c in zip([5, 7, 9], ['0.5','0.6','0.8']):

    # now calculate the observed frame spectra
    
    sed.get_fnu(cosmo, z, igm=Madau96())  # generate observed frame spectra

    ax.plot(sed.obslam/1E4, sed.fnu, lw=1, c=c, zorder=0, alpha = 0.5) # plot SED
   
    fluxes = sed.get_broadband_fluxes(fc)


    # plot fluxes
    for filter, colour in zip(fc, colours):

        f = filter.filter_code

        # ax.scatter(np.log10(filter.pivwv())-4, np.log10(fluxes[f]), c=colour, s=10, zorder=2)
        ax.scatter(filter.pivwv()/1E4, fluxes[f], c=colour, s=10, zorder=2)

        print(filter.pivwv()/1E4, fluxes[f])

    # plot colours (as connecting lines)
    for f1, f2, colour in zip(fc.filter_codes[:-1], fc.filter_codes[1:], colour_colours):

        # x = [np.log10(fc[f1].pivwv())-4, np.log10(fc[f2].pivwv())-4]
        # y = [np.log10(fluxes[f1]), np.log10(fluxes[f2])]
        # ax.plot(x,y, c=colour, lw=1, alpha=0.5, zorder=1)

        x = [fc[f1].pivwv()/1E4, fc[f2].pivwv()/1E4]
        y = [fluxes[f1], fluxes[f2]]
        ax.plot(x,y, c=colour, lw=1, alpha=0.5, zorder=1)


# plot filter transmission function
for filter, colour in zip(fc, colours):
    f = filter.filter_code
    filt_ax.plot(filter.lam/1E4, filter.t, c=colour, lw=2, alpha = 0.5)

    offset = 0.0
    if f == 'JWST/NIRCam.F410M': offset = 0.1

    filt_ax.text(filter.pivwv()/1E4, np.max(filter.t)*1.1 + offset, f.split('.')[-1], c=colour, fontsize=6, ha='center')
        



xlim = [0.9,6.3]

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(xlim)
ax.set_ylim([5, 70])

ax.set_yticks([6,7,8,9,10,20,30,40,50,60])
ax.set_yticklabels([6,7,8,9,10,20,30,40,50,60])



filt_ax.set_xlim(xlim)
filt_ax.set_xscale('log')
filt_ax.set_xticks([1,2,3,4,5,6])
filt_ax.set_xticklabels([1,2,3,4,5,6])
filt_ax.set_ylim([0.01,0.75])
filt_ax.set_yticks([])

# ax.legend(fontsize=8)

filt_ax.set_xlabel(r'$\rm \lambda/\mu m $')
ax.set_ylabel(r'$\rm f_{\nu}/nJy\cdot 10^{8}\ M_{\odot}$')

ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)


fig.savefig(f'figs/sed.pdf')



fig.clf()


