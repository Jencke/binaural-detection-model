"""Calculate the Results for Robinson 1963."""
from ..model import find_snr
import numpy as np
import pandas as pds
import os
from ..model_helpers import band2rms


def calc(rho_hat=0.92, bin_noise=0.31, mon_noise=0.76, calc=True,
                   save=True):
    """Run the experiment."""
    if not calc:
        fname = 'data.h5' if save is True else save
        all_data = pds.read_hdf(fname, 'robinson1963')
        return all_data

    bw = 900

    rel_tone_level = 69.2           # dB
    noise_band_level = 50           # dB
    noise_level = band2rms(noise_band_level, bw)

    path, _ = os.path.split(__file__)
    file = os.path.join(path, 'exp_data.csv')
    exp_data = pds.read_csv(file, header=0, skipinitialspace=True)

    exp_data['Spi'] = rel_tone_level - exp_data['Spi']
    exp_data['S0'] = rel_tone_level - exp_data['S0']
    exp_data['Spi'] -= noise_level
    exp_data['S0'] -= noise_level

    all_data = pds.DataFrame()
    for i, df in exp_data.iterrows():
        spi_dict = {'rho': df.rho,
                    'SNR': df.Spi,
                    'model': False,
                    'exp': 'Spi'}
        s0_dict = {'rho': -df.rho,
                   'SNR': df.S0,
                   'model': False,
                   'exp': 'S0'}
        all_data = all_data.append(pds.DataFrame([spi_dict]))
        all_data = all_data.append(pds.DataFrame([s0_dict]))

    goal_dprime = 1.5

    rho_vals = np.sort(all_data.rho.unique())
    for i_rho, rho in enumerate(rho_vals):
        t_param = dict(noise_bw=bw, n_corr=rho, s_ipd=np.pi)
        r_param = dict(noise_bw=bw, n_corr=rho, s_ipd=np.pi)
        snr_spi = find_snr(bin_noise=bin_noise,
                           rho_hat=rho_hat,
                           mon_noise=mon_noise,
                           goal_dprime=goal_dprime,
                           target_param=t_param,
                           reference_param=r_param)
        t_param = dict(noise_bw=bw, n_corr=rho, s_ipd=0)
        r_param = dict(noise_bw=bw, n_corr=rho, s_ipd=0)
        snr_s0 = find_snr(bin_noise=bin_noise,
                          rho_hat=rho_hat,
                          mon_noise=mon_noise,
                          goal_dprime=goal_dprime,
                          target_param=t_param,
                          reference_param=r_param)
        spi_dict = dict(rho=rho, SNR=10*np.log10(snr_spi[0]), model=True, exp='Spi',
                        bin_noise=bin_noise, mon_noise=mon_noise, rho_hat=rho_hat)
        s0_dict = dict(rho=rho, SNR=10*np.log10(snr_s0[0]), model=True, exp='S0',
                       bin_noise=bin_noise, mon_noise=mon_noise, rho_hat=rho_hat)
        all_data = all_data.append(pds.DataFrame([spi_dict]))
        all_data = all_data.append(pds.DataFrame([s0_dict]))

    if save:
        fname = 'data.h5' if save is True else save
        all_data.to_hdf(fname, 'robinson1963')
    return all_data
