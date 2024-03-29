
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cmasher as cmr
import h5py
from astropy.table import Table
from scipy.stats import binned_statistic

from synthesizer.filters import SVOFilterCollection
from flags.pz import eazy

np.random.seed(42)


def fnu_to_m(fnu):
    """ Convert fnu/nJy to AB magnitude """
    return -2.5*np.log10(fnu/1E9) + 8.9  # -- assumes flux in nJy


def m_to_fnu(m):
    """ Convert AB magnitude to fnu """
    return 1E9 * 10**(-0.4*(m - 8.9))  # -- flux returned nJy


def create_test_observations(N=10):

    from flare_lf.lightcones import Lightcone

    lightcone = Lightcone('jaguar')

    m_ref = lightcone.m[f_ref]

    z = lightcone.z

    s = (z > 4) & (z < 15) & (m_ref < 28.5) & (m_ref > 22)

    with h5py.File(f'data/test.h5', 'w') as hf:

        hf[f'input/z'] = lightcone.z[s][:N]

        for f in filters:

            hf[f'input/{f}/flux'] = lightcone.fnu[f][s][:N]
            hf[f'obs/{f}/flux'] = lightcone.fnu[f][s][:N] + depths_fnu[f]*np.random.normal(size=N)
            hf[f'obs/{f}/flux_err'] = depths_fnu[f] * np.ones(N)


def create_synthetic_observations(model):

    from flare_lf.lightcones import Lightcone

    lightcone = Lightcone(model)

    m_ref = lightcone.m[f_ref]

    z = lightcone.z

    s = (z > 4) & (z < 15) & (m_ref < 28.5) & (m_ref > 22)

    N = np.sum(s)
    print(model, N)

    with h5py.File(f'data/{model}.h5', 'w') as hf:

        hf[f'input/z'] = lightcone.z[s]

        for f in filters:

            hf[f'input/{f}/flux'] = lightcone.fnu[f][s]

            # need to modify this to allow flux dependent noise.

            hf[f'obs/{f}/flux'] = lightcone.fnu[f][s] + depths_fnu[f]*np.random.normal(size=N)
            hf[f'obs/{f}/flux_err'] = depths_fnu[f] * np.ones(N)


def run_eazy(model):

    with h5py.File(f'data/{model}.h5', 'a') as hf:

        id = model

        template_set = 'Larson22'

        # --- initialise EAZY fitter
        pz = eazy.Eazy(id, filter_collection, create_POFZ_FILE=True)

        pz.params['TEMPLATES_FILE'] = f'templates/{template_set}.spectra.param'

        # --- create input catalogue from HDF5 object
        pz.create_input_catalogue_from_HDF5(hf)
        pz.run()

        try:
            del hf['pz/eazy/'+template_set]
        except:
            print('nothing to delete')

        eazy.append_EAZY_output_to_HDF5(
            f'EAZY/outputs/{id}', hf, read_pz=False, read_template_norm=False, get_integrals=True, group_name='pz/eazy/'+template_set)

        # --- need to extract P(z) limits.


if __name__ == "__main__":

    f_ref = 'JWST/NIRCam.F277W'

    # models = ['flares', 'scsam', 'jaguar', 'dream']
    models = ['test', 'scsam', 'jaguar', 'dream', 'flares']

    filters = []
    filters += [f'HST/ACS_WFC.{f}' for f in ['F606W', 'F814W']]
    filters += [f'HST/WFC3_IR.{f}' for f in ['F105W', 'F125W', 'F160W']]
    filters += [f'JWST/NIRCam.{f}' for f in ['F115W',
                                             'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']]

    filter_collection = SVOFilterCollection(filters)  # necessary for EAZY

    depths_m = {  # taken from CEERS Key Paper
        'HST/ACS_WFC.F606W': 28.6,
        'HST/ACS_WFC.F814W': 28.3,
        'HST/WFC3_IR.F105W': 27.1,
        'HST/WFC3_IR.F125W': 27.3,
        'HST/WFC3_IR.F140W': 26.7,
        'HST/WFC3_IR.F160W': 27.4,
        'JWST/NIRCam.F115W': 29.2,
        'JWST/NIRCam.F150W': 29.0,
        'JWST/NIRCam.F200W': 29.2,
        'JWST/NIRCam.F277W': 29.2,
        'JWST/NIRCam.F356W': 29.2,
        'JWST/NIRCam.F410M': 28.4,
        'JWST/NIRCam.F444W': 28.6,
    }

    depths_fnu = {f: m_to_fnu(depths_m[f])/5. for f in filters}  # 1\sigma

    if len(sys.argv) > 1:
        model = models[int(sys.argv[1])-1]
    else:
        model = 'test'

    print(model)
    # create_test_observations()
    # create_synthetic_observations(model)
    run_eazy(model)
