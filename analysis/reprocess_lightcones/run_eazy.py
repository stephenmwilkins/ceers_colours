
import numpy as np
import sys
import h5py
import yaml
from synthesizer.filters import FilterCollection
from flags.pz import eazy


if __name__ == "__main__":

    # list of models to run over
    models = ['scsam1000_basic']

    # define path to eazy
    path_to_eazy = '/Users/sw376/Dropbox/Research/software/eazy-photoz'
    
    # define template set
    template_set = 'Larson22'

    # load eazy parameters
    eazy_params = yaml.safe_load(open('eazy_params.yaml','r'))
    eazy_params['TEMPLATES_FILE'] = f'templates/{template_set}.spectra.param'


    # loop over models
    for model in models:

        # get all filters to make EAZY filter res file

        with h5py.File(f'data/{model}.h5', 'r') as hf:

            # get list of filters
            filters = hf.attrs['filters']

            # initialise filter collection
            filter_collection = FilterCollection(filters, new_lam = np.arange(1000.,100000.,1000.)) # new_lam doesn't do anything here since we use the original wavelength grid

            # initialise EAZY fitter
            pz = eazy.Eazy(model, filter_collection, params = eazy_params, path_to_eazy = path_to_eazy, create_POFZ_FILE=False)

            # create input catalogue from HDF5 object
            pz.create_input_catalogue_from_HDF5(hf, flux_path = lambda f: f'fnu/{f}', flux_err_path = lambda f: f'fnu_err/{f}')
        
        # run eazy
        pz.run()

        with h5py.File(f'data/{model}.h5', 'a') as hf:

            eazy_group = 'pz/eazy/'+template_set
        
            # delete photometric redshifts if they already exist
            if eazy_group in hf:
                del hf[eazy_group]
        
            eazy.append_EAZY_output_to_HDF5(
                f'eazy/outputs/{model}', hf, read_pz=False, read_template_norm=False, get_integrals=False, group_name='pz/eazy/'+template_set)

            # --- need to extract P(z) limits.