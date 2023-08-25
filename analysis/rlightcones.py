



import h5py
import numpy as np
from astropy.table import Table, vstack

from synthesizer.utils import m_to_fnu, fnu_to_m


simulation_data_dir = f'/Users/sw376/Dropbox/Research/data/simulations'




class Lightcone:

    """
    This is for reading reprocessed lightcones
    
    """

    def __init__(self, model, noise_model, data_dir='../reprocess_lightcones/data'):

        if len(model.split('_'))>1:

            model, N = model.split('_')
            filename = f'{data_dir}/{model}{N}_{noise_model}.h5'

        else:

            N = False
            filename = f'{data_dir}/{model}_{noise_model}.h5'




        self.name = model
        self.noise_model = noise_model

        if model == 'flares':
            self.area = 1. # arcmin2
        if model == 'scsam':
            self.area = 782. # arcmin2
        if model == 'jaguar':
            self.area = 11.*11. # arcmin2
        if model == 'dream':
            self.area = 3600. # arcmin2

       

        if model == 'flares':
            self.name = 'FLARES'
        if model == 'scsam':
            self.name = 'SCSAM'
        if model == 'jaguar':
            self.name = 'JAGUAR'
        if model == 'dream':
            self.name = 'DREAM'

        with h5py.File(filename, 'r') as hf:

            # if sub-sampled reduce area
            if N:
                self.area *= hf.attrs['fraction']

            self.filters = hf.attrs['filters']
            print(self.filters)

            self.m = {}
            self.fnu = {}
            self.fnu_err = {}
            for filter in self.filters: 
                self.fnu[filter] = hf[f'fnu/{filter}'][()]
                self.fnu_err[filter] = hf[f'fnu_err/{filter}'][()]
                self.m[filter] = fnu_to_m(self.fnu[filter])

            self.z = hf['z'][()] # input redshift
            self.N = len(self.z)

            # if no noise the photo-z is the true z
            if noise_model == 'no':
                self.pz = self.z
            else:
                # add pz
                hf.visit(print)
                self.pz = hf['pz/eazy/z_a'][()]
       
    def get_selection(self):

        # exclude sources with photometric redshifts of 0.0 as that probably means eazy failed
        s = self.pz != 0.0

        # flux cut
        s = s & (self.fnu['JWST/NIRCam.F277W']>50.)

        

        return s
       


    def Nz(self, area, m_limit, detection_filter, redshift_limits = [5, 15]):

        flux_limit = m_to_fnu(m_limit)
        bin_edges = np.arange(*redshift_limits, 0.01)
        bin_centres = 0.5*(bin_edges[:-1]+bin_edges[1:])

        m = self.m[detection_filter] # magnitudes
        s = (m<m_limit)&(self.z>redshift_limits[0])&(self.z<redshift_limits[1])

        N, _  = np.histogram(self.z[s], bins = bin_edges)
        N = N*area/self.area

        return bin_centres, N


    def CNz(self, area, m_limit, detection_filter, redshift_limits = [5, 15]):

        bin_centres, N = self.Nz(area, m_limit, detection_filter, redshift_limits = redshift_limits)

        CN = np.cumsum(N[::-1])
        return bin_centres[::-1], CN


# if __name__ == "__main__":
#
#     lc = Lightcone('dream')
