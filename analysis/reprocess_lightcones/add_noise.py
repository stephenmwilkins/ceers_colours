



import numpy as np
import h5py
from synthesizer.utils import m_to_fnu, fnu_to_m

# add directory above to Python Path 
import sys
sys.path.append('../')

import lightcones


if __name__ == "__main__":


    noise_model = 'no' # no noise 
    # noise_model = 'flat' # apply flat noise (i.e. same in every band)
    # noise_model = 'basic' # apply uniform noise using stated depths of the images 
    noise_model = 'linear' # apply noise as a linear function of the log10(flux)

    models = ['flares','jaguar','dream', 'scsam']
    models = ['scsam','jaguar','dream']
    models = ['jaguar','dream']

    for model in models:

        # read lightcone
        lightcone = lightcones.Lightcone(model)

        # print the number of galaxies in the lightcone
        print('Total number available:', lightcone.N)

        # limit to apply

        limit_filter, limit_flux = 'JWST/NIRCam.F277W', 25 # nJy

        s = lightcone.fnu[limit_filter] > limit_flux

        N = np.sum(s)

        print('Number selectable:', N)

        # noise_model

        

        

        # used for both basic and linear. For linear this is to provide HST noise
        # taken from CEERS Key Paper, 5sigma depths
        depths_m5 = {  
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

        depths_fnu = {f: m_to_fnu(depths_m5[f]).to('nJy').value/5. for f in lightcone.filters}  # 1\sigma


        if noise_model == 'linear':

            linear_fit = {'JWST/NIRCam.F115W':[0.2675434889093099,0.1190720404415787], 'JWST/NIRCam.F150W':[0.27666798324651953,0.21630720958354066],'JWST/NIRCam.F200W':[0.27835152031590255,0.13058300028405245],'JWST/NIRCam.F277W':[0.294340121054273,-0.03548792067026589],'JWST/NIRCam.F356W':[0.2564473395677601,-0.03418413637403682],'JWST/NIRCam.F410M':[0.21963887365346557,0.3624582113229585],'JWST/NIRCam.F444W':[0.21095758346624174,0.277207138454342]}


        filename = f'{model}_{noise_model}'

        with h5py.File(f'data/{filename}.h5', 'w') as hf:



            # save the input redshift
            hf[f'z'] = lightcone.z[s]
            hf.attrs['filters'] = lightcone.filters

            for f in lightcone.filters:

                # no noise
                if noise_model == 'no':

                    # just use intrinsic
                    hf[f'fnu/{f}'] = lightcone.fnu[f][s]

                    hf[f'fnu_err/{f}'] = np.zeros(N)


                if noise_model == 'basic':

                    # add gaussian noise to the lightcone
                    hf[f'fnu/{f}'] = lightcone.fnu[f][s] + depths_fnu[f]*np.random.normal(size=N)

                    # define error 
                    hf[f'fnu_err/{f}'] = depths_fnu[f] * np.ones(N)

                if noise_model == 'linear':

                    if f.split('/')[0] == 'JWST':

                        

                        b, c = linear_fit[f]

                        sigma = 10**(b*np.log10(lightcone.fnu[f][s])+c) 

                        # add gaussian noise to the lightcone
                        hf[f'fnu/{f}'] = lightcone.fnu[f][s] + sigma*np.random.normal(size=N)

                        # define error 
                        hf[f'fnu_err/{f}'] = sigma

                        print(f, 'JWST', f'{min(sigma):.2f} {max(sigma):.2f}')

                    else:
                        
                        print(f, 'HST')

                        # add gaussian noise to the lightcone
                        hf[f'fnu/{f}'] = lightcone.fnu[f][s] + depths_fnu[f]*np.random.normal(size=N)

                        # define error 
                        hf[f'fnu_err/{f}'] = depths_fnu[f] * np.ones(N)


            