
import numpy as np

from astropy.io import fits
from astropy.table import Table, vstack

import matplotlib as mpl
import matplotlib.pyplot as plt
import cmasher as cmr

from flare.photom import flux_to_m as fnu_to_m

import flare.plt as fplt



class ReadCEERS:

    def __init__(self, reference_filter = 'F200', version = '0.2'):

        self.version = version
        self.area = 34. # arcmin2

        data_dir = '/Users/stephenwilkins/Dropbox/Research/data/images/jwst/ceers/cats'

        for i, field in enumerate([1,2,3,6]):

            cat_name = f'CEERS_NIRCam{field}_v{version}'
            filename = f'{data_dir}/{cat_name}'

            if i == 0:
                photo_cat = Table.read(f'{filename}_photom.fits')
                zphot_cat = Table.read(f'{filename}_zphot.fits')

                photo_cat['ZA'] = zphot_cat['ZA'][0]


            else:

                photo_cat_ = Table.read(f'{filename}_photom.fits')
                zphot_cat_ = Table.read(f'{filename}_zphot.fits')

                photo_cat_['ZA'] = zphot_cat_['ZA'][0]

                photo_cat = vstack([photo_cat, photo_cat_])
                # zphot_cat = vstack([zphot_cat, zphot_cat_])

        self.z = photo_cat['ZA']

        self.photo_cat = photo_cat

        self.m = {}
        self.flux = {}
        self.flux_err = {}

        for f in ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']:
            self.flux[f] = photo_cat[f[:-1]]
            self.flux_err[f] = photo_cat['D'+f[:-1]]
            self.m[f] = fnu_to_m(photo_cat[f[:-1]])


    def get_selection(self, reference_filter = 'F200W', sn_cut = 10, magnitude_limits = (24, 28), redshift_limits = (4, 10)):

        sn = self.flux[reference_filter]/self.flux_err[reference_filter]

        s = (sn>sn_cut)&(self.z>redshift_limits[0])&(self.z<redshift_limits[1])&(self.m[reference_filter]>magnitude_limits[0])&(self.m[reference_filter]<magnitude_limits[1])

        print(f'number of sources selected: {np.sum(s)}')

        return s
