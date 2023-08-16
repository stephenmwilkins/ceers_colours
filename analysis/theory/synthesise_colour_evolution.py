
import pickle
import numpy as np
from synthesizer.filters import FilterCollection
from synthesizer.grid import Grid
from synthesizer.parametric.sfzh import SFH, ZH, generate_sfzh
from synthesizer.parametric.galaxy import Galaxy
from unyt import yr, Myr
from synthesizer.igm import Madau96, Inoue14
from astropy.cosmology import Planck18 as cosmo


# define choise of SPS model and initial mass function (IMF)

grid_name = "bpass-2.2.1-bin_chabrier03-0.1,300.0_cloudy-c17.03"
grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids'
grid = Grid(grid_name, grid_dir=grid_dir)

# define filters
filter_codes = [f'JWST/NIRCam.{f}' for f in ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']]  # define a list of filter codes
fc = FilterCollection(filter_codes, new_lam=grid.lam)
colours = [[f1, f2] for f1, f2 in zip(filter_codes[:-1], filter_codes[1:])]


for duration in [10, 100]:

    # define the parameters of the star formation and metal enrichment histories
    sfh_p = {'duration': duration * Myr}
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

    stellar_sed = galaxy.get_spectra_stellar(grid)
    intrinsic_sed = galaxy.get_spectra_intrinsic(grid)


    # define array of redshifts
    redshifts = np.arange(4., 10.,0.01)

    # create outputs
    output = {}
    output['z'] = redshifts
    for t in ['stellar','intrinsic']:
        output[t] = {}
        for f in filter_codes:
            output[t][f] = np.empty(len(redshifts))

    # loop over redshifts
    for i, z in enumerate(redshifts):
        
        print(f'{z:.2f}')

        # generate observed frame **stellar** spectra, fluxes, and save
        stellar_sed.get_fnu(cosmo, z, igm=Madau96())  
        fluxes = stellar_sed.get_broadband_fluxes(fc)
        for f in filter_codes:
            output['stellar'][f][i] = fluxes[f].value

        # generate observed frame **intrinsic** spectra, fluxes, and save
        intrinsic_sed.get_fnu(cosmo, z, igm=Madau96())  
        fluxes = intrinsic_sed.get_broadband_fluxes(fc)
        for f in filter_codes:
            output['intrinsic'][f][i] = fluxes[f].value


    pickle.dump(output, open(f'data/{duration}Myr.pck','wb'))