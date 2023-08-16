
import h5py
import numpy as np
from astropy.table import Table, vstack

from synthesizer.utils import m_to_fnu, fnu_to_m


simulation_data_dir = f'/Users/sw376/Dropbox/Research/data/simulations'


models = ['scsam']

# filters = ['Webb.NIRCAM.F277W']
# filters = [f'Webb.NIRCAM.{filter}' for filter in ['F115W','F150W', 'F200W','F277W','F356W','F410M','F444W']]

filters = []
filters += [f'HST/ACS_WFC.{f}' for f in ['F606W', 'F814W']]
filters += [f'HST/WFC3_IR.{f}' for f in ['F105W', 'F125W', 'F160W']]
filters += [f'JWST/NIRCam.{f}' for f in ['F115W','F150W', 'F200W','F277W','F356W','F410M','F444W']]



class Lightcone:

    def __init__(self, model):


        print('-'*30)
        print(model)

        if model == 'flares':

            self.name = r'FLARES'
            self.area = 1. # arcmin2
            self.area_deg = self.area/3600.

            self.data = h5py.File(f'{simulation_data_dir}/flares/colours/DustModelI.h5','r')

            filter_dict = {}
            # filter_dict = filter_dict | {f'HST/WFC3_IR.{f}': f'HST.WCF3.{f}' for f in ['F105W', 'F125W', 'F160W']}
            filter_dict = filter_dict | {f'JWST/NIRCam.{f}': f'Webb.NIRCAM.{f}' for f in ['F115W', 'F150W', 'F200W','F277W','F356W','F410M','F444W']}


            self.m = {}
            self.fnu = {}
            for filter in filters:
                if filter in filter_dict.keys():
                    f = filter_dict[filter]
                    self.fnu[filter] = self.data[f][()]
                    self.m[filter] = fnu_to_m(self.fnu[filter])

            self.z = self.data['z'][()]


        elif model == 'scsam':
            filename = 'JWST_phot_dust'
            simulation = 'SCSAM/egs.0'
            self.name = r'Santa\ Cruz\ SAM'
            self.area = 782. # arcmin2
            self.area_deg = self.area/3600.
            self.data = h5py.File(f'{simulation_data_dir}/{simulation}/{filename}.h5', 'r')

            self.m = {}
            self.fnu = {}

            filter_dict = {}
            filter_dict = filter_dict | {f'HST/ACS_WFC.{f}': f'acs{f.lower()}_dust' for f in ['F606W', 'F814W']}
            filter_dict = filter_dict | {f'HST/WFC3_IR.{f}': f'wfc3{f.lower()}_dust' for f in ['F105W', 'F125W', 'F160W']}
            filter_dict = filter_dict | {f'JWST/NIRCam.{f}': f'NIRCam_{f}_dust' for f in ['F115W', 'F150W', 'F200W','F277W','F356W','F410M','F444W']}

            for filter in filters:
                if filter in filter_dict.keys():
                    f = filter_dict[filter]
                    self.m[filter] = self.data[f][()]
                    self.fnu[filter] = m_to_fnu(self.m[filter])

            self.z = self.data['redshift'][()]


        elif model == 'jaguar':
            self.name = r'JAGUAR'
            self.area = 11.*11. # arcmin2
            data_SF = Table.read(f'{simulation_data_dir}/JAGUAR/JADES_SF_mock_r1_v1.2.fits')
            data_Q = Table.read(f'{simulation_data_dir}/JAGUAR/JADES_Q_mock_r1_v1.2.fits')

            self.data = vstack([data_SF, data_Q])

            # print(self.data.keys())

            filter_dict = {}
            filter_dict = filter_dict | {f'HST/ACS_WFC.{f}': f'HST_{f}_fnu' for f in ['F606W', 'F814W']}
            filter_dict = filter_dict | {f'HST/WFC3_IR.{f}': f'HST_{f}_fnu' for f in ['F105W', 'F125W', 'F160W']}
            filter_dict = filter_dict | {f'JWST/NIRCam.{f}': f'NRC_{f}_fnu' for f in ['F115W', 'F150W', 'F200W','F277W','F356W','F410M','F444W']}

            self.m = {}
            self.fnu = {}
            for filter in filters:
                f = filter_dict[filter]
                self.fnu[filter] = self.data[f].data
                self.m[filter] = -2.5*np.log10(self.fnu[filter]/1E9) + 8.9
            self.z = self.data['redshift'].data

        elif model == 'dream':

            self.name = r'DREaM'
            self.area = 3600. # arcmin2

            self.data = h5py.File(f'{simulation_data_dir}/DREaM/DREaM_photo.h5', 'r')

            filter_dict = {}
            filter_dict = filter_dict | {f'HST/ACS_WFC.{f}': f'WFC_ACS_{f}' for f in ['F606W', 'F814W']}
            filter_dict = filter_dict | {f'HST/WFC3_IR.{f}': f'WFC3_IR_{f}' for f in ['F105W', 'F125W', 'F160W']}
            filter_dict = filter_dict | {f'JWST/NIRCam.{f}': f'JWST_{f}' for f in ['F115W', 'F150W', 'F200W','F277W','F356W','F410M','F444W']}

            self.m = {}
            self.fnu = {}
            for filter in filters:
                f = filter_dict[filter]
                self.m[filter] = self.data[f][()]
                self.fnu[filter] = m_to_fnu(self.m[filter])

            self.z = self.data['redshift'][()]


        else:

            print('ERROR: model name not recognised')

        print('initialised model')


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
