



import numpy as np
import h5py
from synthesizer.utils import m_to_fnu, fnu_to_m

# add directory above to Python Path 
import sys
sys.path.append('../')




if __name__ == "__main__":

    N = 100

    
    noise_model = 'basic' # apply uniform noise using stated depths of the images 
    # noise_model = 'linear' # apply noise as a linear function of the log10(flux)

    models = ['flares','jaguar','dream', 'scsam']
    models = ['scsam','jaguar','dream']
    # models = ['jaguar']

    for model in models:

        print(model)
        filename = f'{model}_{noise_model}'

        # open resampled 
        with h5py.File(f'data/{model}{N}_{noise_model}.h5', 'w') as hf:

            with h5py.File(f'data/{filename}.h5', 'r') as hfin:



                Ntot = len(hfin[f'z'])
                fraction = N/Ntot

                ids = np.arange(Ntot)

                print(f'fraction sampled: {fraction:.3f}')
                hf.attrs['fraction'] = fraction
                hf.attrs['filters'] = hfin.attrs['filters']
                s = np.sort(np.random.choice(ids,size=N,replace=False))

                # save the input redshift
                hf[f'z'] = hfin[f'z'][s]

                for f in hfin.attrs['filters']:

                    # just use intrinsic
                    hf[f'fnu/{f}'] = hfin[f'fnu/{f}'][s]
                    hf[f'fnu_err/{f}'] = hfin[f'fnu_err/{f}'][s]

                