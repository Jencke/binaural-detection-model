"""Calculate model results for van de Par 1999."""
import numpy as np
import os
from ..model import find_snr
import pandas as pds



def calc(rho_hat=0.97, bin_noise=0.38, mon_noise=0.76, calc=True,
         save=True):
    """Run the Experiment."""
    if not calc:
        fname = 'data.h5' if save is True else save
        all_data = pds.read_hdf(fname, 'vandepar1999')
        return all_data

    path, _ = os.path.split(__file__)
    file = os.path.join(path, 'exp_data.csv')
    exp_data = pds.read_csv(file, header=0, skipinitialspace=True)

    all_data = pds.DataFrame()
    for i, df in exp_data.iterrows():
        n0spi_dict = {'bw': df.bw,
                      'SNR': df.N0Spi,
                      'model': False,
                      'exp': 'N0Spi'}
        n0s0_dict = {'bw': df.bw,
                     'SNR': df.N0S0,
                     'model': False,
                     'exp': 'N0S0'}
        npis0_dict = {'bw': df.bw,
                      'SNR': df.NpiS0,
                      'model': False,
                      'exp': 'NpiS0'}
        all_data = all_data.append(pds.DataFrame([n0spi_dict]))
        all_data = all_data.append(pds.DataFrame([n0s0_dict]))
        all_data = all_data.append(pds.DataFrame([npis0_dict]))

    bw_values = [5, 10, 25, 50, 100, 250, 500, 1000]
    goal_dprime = 0.77

    for i_bw, bw in enumerate(bw_values):
        t_param = dict(noise_bw=bw,
                       s_ipd=np.pi)
        r_param = dict(noise_bw=bw,
                       s_ipd=np.pi)
        snr_n0spi = find_snr(bin_noise=bin_noise,
                             rho_hat=rho_hat,
                             mon_noise=mon_noise,
                             goal_dprime=goal_dprime,
                             target_param=t_param,
                             reference_param=r_param,
                             ftol=1e-6
                             )
        t_param = dict(noise_bw=bw,
                       s_ipd=0)
        r_param = dict(noise_bw=bw,
                       s_ipd=0)
        snr_n0s0 = find_snr(bin_noise=bin_noise,
                            rho_hat=rho_hat,
                            mon_noise=mon_noise,
                            goal_dprime=goal_dprime,
                            target_param=t_param,
                            reference_param=r_param,
                            ftol=1e-6
                            )
        n0spi_dict = dict(SNR=10*np.log10(snr_n0spi[0]),
                          model=True, exp='N0Spi', bw=bw,
                          bin_noise=bin_noise,
                          mon_noise=mon_noise,
                          rho_hat=rho_hat)
        npis0_dict = dict(SNR=10*np.log10(snr_n0spi[0]),
                          model=True, exp='NpiS0', bw=bw,
                          bin_noise=bin_noise,
                          mon_noise=mon_noise,
                          rho_hat=rho_hat)
        n0s0_dict = dict(SNR=10*np.log10(snr_n0s0[0]), model=True,
                         exp='N0S0', bw=bw, bin_noise=bin_noise,
                         mon_noise=mon_noise, rho_hat=rho_hat)
        all_data = all_data.append(pds.DataFrame([n0spi_dict]))
        all_data = all_data.append(pds.DataFrame([npis0_dict]))
        all_data = all_data.append(pds.DataFrame([n0s0_dict]))

    if save:
        fname = 'data.h5' if save is True else save
        all_data.to_hdf(fname, 'vandepar1999')
    return all_data
#
