
import numpy as np
import sys
import h5py
import yaml
from synthesizer.filters import FilterCollection
from flags.pz import eazy
import sys


if __name__ == "__main__":

    jobs = 99 # actually need to run jobs + 1

    model = sys.argv[1] # model name (e.g. scsam_basic)
    i = int(sys.argv[2])-1 # iteration

    # define path to eazy
    path_to_eazy = '/Users/sw376/Dropbox/Research/software/eazy-photoz'
    
    # define template set
    template_set = 'Larson22'

    # load eazy parameters
    eazy_params = yaml.safe_load(open('eazy_params.yaml','r'))
    eazy_params['TEMPLATES_FILE'] = f'templates/{template_set}.spectra.param'

    

    with h5py.File(f'data/{model}.h5', 'r') as hf:

        # get list of filters
        filters = hf.attrs['filters']

        Ntot = len(hf[f'fnu/{filters[0]}'])

        N = int(Ntot/jobs)

        start = N*i
        end = N*(i+1)

        if end > Ntot: end = Ntot

        print(i, start, end, Ntot)

        flux_dict = {}
        err_dict = {}

        for f in filters:
            flux_dict[f] = hf[f'fnu/{f}'][start:end]
            err_dict[f] = hf[f'fnu_err/{f}'][start:end]


    # initialise filter collection
    filter_collection = FilterCollection(filters, new_lam = np.arange(1000.,100000.,1000.)) # new_lam doesn't do anything here since we use the original wavelength grid

    # initialise EAZY fitter
    pz = eazy.Eazy(model, filter_collection, params = eazy_params, path_to_eazy = path_to_eazy, create_POFZ_FILE=False)

    # create input catalogue from HDF5 object
    pz.create_input_catalogue_from_dict(flux_dict, err_dict)

    # run eazy
    pz.run()

    # need to be re-combined later