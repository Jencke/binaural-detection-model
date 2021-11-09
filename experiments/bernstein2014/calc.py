"""Calculate model results for Bernstein 2014."""

import numpy as np
from ..model import find_snr
import pandas as pds
import os

def calc(rho_hat=0.97, bin_noise=0.54, mon_noise=0.76, calc=True,
        save=True):
    """Run the experiment."""

    if not calc:
        fname = 'data.h5' if save is True else save
        all_data = pds.read_hdf(fname, 'bernstein2014')
        return all_data

    path, _ = os.path.split(__file__)
    file = os.path.join(path, 'exp_data.csv')
    exp_data = pds.read_csv(file, header=0, skipinitialspace=True)

    bw_values = [25, 100, 300, 900]
    rho_values = exp_data.rho
    rho_values = np.unique(np.sort(rho_values))

    goal_dprime = 0.770

    thr = np.zeros([len(bw_values), len(rho_values)])
    for i_bw, bw in enumerate(bw_values):
        for i_rho, rho in enumerate(rho_values):
            t_param = dict(noise_bw=bw,
                           n_corr=rho,
                           s_ipd=np.pi)
            r_param = dict(noise_bw=bw,
                           n_corr=rho,
                           s_ipd=np.pi)
            snr_n0spi = find_snr(bin_noise=bin_noise,
                                 rho_hat=rho_hat,
                                 mon_noise=mon_noise,
                                 goal_dprime=goal_dprime,
                                 target_param=t_param,
                                 reference_param=r_param
                                 )
            thr[i_bw, i_rho] = 10 * np.log10(snr_n0spi)

    # Combine experimental with model data
    all_data = pds.DataFrame()
    for i_bw, bw in enumerate(bw_values):
        all_data = all_data.append(pds.DataFrame({'rho': rho_values,
                                                  'SNR': thr[i_bw],
                                                  'bw': bw,
                                                  'model': True,
                                                  'exp': 'n_rho_s_pi',
                                                  'rho_hat': rho_hat,
                                                  'mon_noise': mon_noise,
                                                  'bin_noise': bin_noise
                                                  }))
        all_data = all_data.append(pds.DataFrame({'rho': exp_data.rho,
                                                  'SNR': exp_data[f"{bw}hz"],
                                                  'bw': bw,
                                                  'model': False,
                                                  'exp': 'n_rho_s_pi'}))

    all_data['paper'] = 'bernstein2014'
    if save:
        fname = 'data.h5' if save is True else save
        all_data.to_hdf(fname, 'bernstein2014')
    return all_data
