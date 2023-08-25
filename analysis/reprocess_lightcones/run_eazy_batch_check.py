
import numpy as np
import sys
import h5py
from astropy.io import ascii


if __name__ == "__main__":

    jobs = 99 # actually need to run jobs + 1

    model = sys.argv[1] # model name (e.g. scsam_basic)
    
    with h5py.File(f'data/{model}.h5', 'a') as hf:

        # get list of filters
        filters = hf.attrs['filters']

        Ntot = len(hf[f'fnu/{filters[0]}'])

        N = int(Ntot/jobs)

        for i in range(jobs):

            with open(f'eazy/outputs/{model}_{i}.zout','r') as f:
                print(i, len(f.readlines()))

           
            




