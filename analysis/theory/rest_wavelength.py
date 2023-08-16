

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from synthesizer.filters import FilterCollection
from synthesizer.grid import Grid
from synthesizer.parametric.sfzh import SFH, ZH, generate_sfzh
from synthesizer.parametric.galaxy import Galaxy
from unyt import yr, Myr



plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')




# predict EW for all the lines

grid_name = "bpass-2.2.1-bin_chabrier03-0.1,300.0_cloudy-c17.03"
grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids'
grid = Grid(grid_name, grid_dir=grid_dir, read_spectra=False, read_lines=True)

# define the parameters of the star formation and metal enrichment histories
sfh_p = {'duration': 10 * Myr}
Z_p = {'log10Z': -2.5}  # can also use linear metallicity e.g. {'Z': 0.01}

# define the functional form of the star formation and metal enrichment histories
sfh = SFH.Constant(sfh_p)  # constant star formation
Zh = ZH.deltaConstant(Z_p)  # constant metallicity

# get the 2D star formation and metal enrichment history for the given SPS grid and print summary.
sfzh = generate_sfzh(grid.log10age, grid.metallicity, sfh, Zh)

# create the Galaxy object and print a summary
galaxy = Galaxy(sfzh)


# define list of lines that we're interested in. Note that we can provide multiples which are automatically summed
line_ids = ['H 1 4862.69A', 'O 3 4958.91A', 'O 3 5006.84A', ['O 3 4958.91A', 'O 3 5006.84A']]

# create the Lines dictionary which contains line objects
lines = galaxy.get_line_intrinsic(grid, grid.line_list)


z = np.arange(5,10,0.01)


# ----- second plot

fig = plt.figure(figsize = (3.5, 4))

left  = 0.15
height = 0.85
bottom = 0.1
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


# --- determine line contamination
for line_name, line in lines.items():

    c = 'k'
    # set lw and alpha based on line EW
    alpha = np.sqrt(line.ew.value/3000)
    lw = np.sqrt(3*line.ew.value/3000)
    print(line_name, line.ew.value)
    ax.axhline(line.wavelength.value, c=c, alpha = alpha, lw=lw, zorder=2)

# --- plot break boundary


ax.axhline(912, c='k', lw=1, ls='--', alpha=1.0, zorder=2)
ax.axhline(1216, c='k', lw=1, ls='--', alpha=1.0, zorder=2)
ax.axhline(3700, c='k', lw=1, ls='--', alpha=1.0, zorder=2)


# define filters
filter_codes = [f'JWST/NIRCam.{f}' for f in ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']]  # define a list of filter codes
fc = FilterCollection(filter_codes)

# filter colours
colours = cm.rainbow(np.linspace(0, 1, len(fc.filters)))

for i,(filter,c) in enumerate(zip(fc, colours)):
    mn, mx = filter.mnmx()
    piv = filter.pivwv()

    d = mn/(1+z)
    u = mx/(1+z)

    ax.fill_between(z, d, u, color=c, alpha = 0.5, label = rf'$\rm {filter.filter_code.split(".")[-1]}$', zorder=0)
    ax.plot(z, d,color=c, lw=1, zorder=1)
    ax.plot(z, u,color=c, lw=1, zorder=1)
    # ax.text(12.1, piv/13., f.split('.')[-1], fontsize=5, c=c)


ax.legend(fontsize=7, loc = 'upper right', labelspacing=0.1)


# ax.set_xlim([5., 12])
ax.set_ylim([100, 9500])
ax.set_xlim([5., 10])
ax.set_xlabel(r'$\rm z$')
ax.set_ylabel(r'$\rm \lambda_{rest}/{\regular \AA}$')

fig.savefig(f'figs/rest_wavelength.pdf')
fig.clf()
