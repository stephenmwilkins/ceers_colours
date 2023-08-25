
import numpy as np
import sys
import h5py



if __name__ == "__main__":

    jobs = 99 # actually need to run jobs + 1

    model = sys.argv[1] # model name (e.g. scsam_basic)
    
    with h5py.File(f'data/{model}.h5', 'a') as hf:

        # get list of filters
        filters = hf.attrs['filters']

        Ntot = len(hf[f'fnu/{filters[0]}'])

        N = int(Ntot/jobs)

        # set up output datasets
        eazy_datasets = ['z_a','z_m1','chi_a','l68','u68','l95','u95','l99','u99','q_z','z_peak','peak_prob','z_mc']

        if 'pz/eazy/' in hf:
            del hf['pz/eazy/']

        for ds in eazy_datasets:
            hf[f'pz/eazy/{ds}'] = np.empty(Ntot)

        for i in range(jobs):

            # open zout file
            zout = ascii.read(f'eazy/output/{model}_{i}.zout')

            start = N*i
            end = N*(i+1)
            if end > Ntot: end = Ntot

            for ds in eazy_datasets:
                hf[f'pz/eazy/{ds}'][start:end] = zout[ds]
            




