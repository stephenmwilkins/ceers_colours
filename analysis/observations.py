

import h5py
import numpy as np
from astropy.cosmology import Planck18 as cosmo
from unyt import unyt_array

from synthesizer.filters import FilterCollection
from synthesizer.utils import flux_to_luminosity, Lnu_to_M, fnu_to_m


# def fnu_to_Lnu(fnu, z, cosmo = Planck18):

#     """
#     Convert flux (in fnu) to luminosity.
#     """

#     # convert to nJy is a unyt array
#     if isinstance(fnu, unyt_array):
#         f = fnu.to('nJy').value
#     # else assume it is in nJy
#     else:
#         f = fnu 
    
#     # convert to W/cm^2/Hz
#     f /= (1E9 * 1E23) 

#     # luminosity distance
#     luminosity_distance = cosmo.luminosity_distance(z).to('cm').value

#     return f*(4.*np.pi*luminosity_distance**2)/(1.+z)


# def Lnu_to_M(Lnu):






class CEERS:

    def __init__(self, data_dir = '../../data', version = 'v3'):

        self.pointings = [1,2,3,4,5,6,7,8,9,10]
        self.area = 90. # arcmin2

        self.cat = {}
        for pointing in self.pointings:
            filename = f'CEERS_NIRCam{pointing}_v0.51.3-CEERS_colours_{version}.h5'
            self.cat[pointing] = h5py.File(f'{data_dir}/{filename}', 'r')
            
        filter_codes = [f'JWST/NIRCam.{f}' for f in ['F115W','F150W', 'F200W','F277W','F356W','F410M','F444W']]
        self.filters = FilterCollection(filter_codes = filter_codes, new_lam = np.arange(1,50000,1))

        # define output arrays
        self.z = np.array([])
        self.id = np.array([])

        self.s = {}
        self.p = np.array([])

        self.fnu = {}
        for filter_code in filter_codes:
            self.fnu[filter_code] = np.array([])

        for pointing in self.pointings:
            
            n = len(self.cat[pointing]['close'][:])

            s = (self.cat[pointing]['close'][:] == False) & (self.cat[pointing]['spurious'][:] == False)

            self.s[pointing] = s

            # get redshift
            z_ = self.cat[pointing]['pz/ceers/ZA'][s]

            # stack
            self.z = np.hstack((self.z,z_))

            p_ = pointing*np.ones(n)
            self.p = np.hstack((self.p,p_[s]))

            # id
            # id = self.cat[pointing][''][s]
            # self.id = np.hstack((self.id,id))


            for filter_code in filter_codes:
                filter_code_ = filter_code.split('.')[-1][1:-1]
                self.fnu[filter_code] = np.hstack((self.fnu[filter_code],self.cat[pointing][f'photom/FLUX_{filter_code_}'][s]))

        self.N = len(self.z)

        # get apparent magnitude
        self.m = {}
        for filter_code in filter_codes:
            self.m[filter_code] = fnu_to_m(self.fnu[filter_code])



    def calculate_M(self, lam=5500.):

        self.Lnu = np.zeros(self.N)
        self.M = np.zeros(self.N)

        for i,z in enumerate(self.z):

            # filter to use to calculate M_1500
            filter_code = self.filters.find_filter(5500., redshift=z).filter_code
            
            fnu = self.fnu[filter_code][i]

            Lnu = flux_to_luminosity(fnu, redshift=z, cosmo=cosmo)
            M = Lnu_to_M(Lnu)

            print(f'{z:.2f} {filter_code} {np.log10(fnu):.2f}, {np.log10(Lnu):.2f} {M:.2f}')

            self.Lnu[i] = Lnu
            self.M[i] = M


            
    def c(self, f1, f2):

        """
        Return an array for a particular colour using SVO filter codes
        """
    
        # get colour
        c = 2.5*np.log10(self.fnu[f2]/self.fnu[f1])

        return c

