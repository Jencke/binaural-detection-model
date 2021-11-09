"""Calculate the Results for Pollack 1959."""
import numpy as np
from ..model import find_rho_jnd
import os
import pandas as pds


def calc(rho_hat=0.92, bin_noise=0.42, mon_noise=0.7, calc=True,
                   save=True, opt=False):
    """Run the model."""
    if not calc:
        fname = 'data.h5' if save is True else save
        all_data = pds.read_hdf(fname, 'pollack1959')
        return all_data

    path, _ = os.path.split(__file__)
    file = os.path.join(path, 'exp_data.csv')
    exp_data = pds.read_csv(file, header=0, skipinitialspace=True)

    all_data = exp_data.copy()
    all_data['model'] = False

    dprimes = [0.36, 0.95, 1.8]
    if not opt:
        ref_corr = 1 - np.logspace(-1.7, 0, 21)
    else:
        ref_corr = []
    ref_corr = np.unique(np.concatenate((ref_corr, all_data['corr'])))
    ref_corr.sort()
    res = np.zeros([len(ref_corr), len(dprimes)])
    for i_d, d in enumerate(dprimes):
        for i_corr, corr in enumerate(ref_corr):
            t_param = dict(noise_bw=1000, snr=0, n_corr=corr)
            r_param = dict(noise_bw=1000, snr=0, n_corr=corr)
            res = find_rho_jnd(bin_noise=bin_noise,
                               rho_hat=rho_hat,
                               mon_noise=mon_noise,
                               target_param=t_param,
                               reference_param=r_param,
                               goal_dprime=d)
            df = pds.DataFrame(dict(jnd=res, corr=corr,
                               dprime=d, model=True,
                               rho_hat=rho_hat,
                               bin_noise=bin_noise,
                               mon_noise=mon_noise))
            all_data = all_data.append(df)
    all_data['paper'] = 'pollack1959'
    if save:
        fname = 'data.h5' if save is True else save
        all_data.to_hdf(fname, 'pollack1959')
    return all_data
