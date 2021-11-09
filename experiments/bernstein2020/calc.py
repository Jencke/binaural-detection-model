"""Calculate the model results from Bernstein 2020."""
import numpy as np
import pandas as pds
import os

from ..model import find_snr


def calc(rho_hat=0.89, bin_noise=0.52, mon_noise=0.93, calc=True,
         save=True, opt=False):
    """Run the experiment."""
    if not calc:
        fname = 'data.h5' if save is True else save
        all_data = pds.read_hdf(fname, 'bernstein2020')
        return all_data

    path, _ = os.path.split(__file__)
    file = os.path.join(path, 'exp_data.csv')
    exp_data = pds.read_csv(file, header=0, skipinitialspace=True)

    data_100 = exp_data[exp_data.bw == 100]
    data_900 = exp_data[exp_data.bw == 900]
    rho_values = [1.0, 0.992, 0.874, 0.498]
    data = pds.DataFrame()
    for i_rho, rho in enumerate(rho_values):
        a = pds.DataFrame({'tau': data_100['tau'],
                           'rho': rho,
                           'SNR': data_100[f'{rho}'].values,
                           'bw': 100,
                           'exp': 'n0srho'})
        b = pds.DataFrame({'tau': data_900['tau'],
                           'rho': rho,
                           'SNR': data_900[f'{rho}'].values,
                           'bw': 900,
                           'exp': 'n0srho'})
        data = data.append([a, b])
    a = pds.DataFrame({'tau': 0,
                       'rho': 1,
                       'SNR': -3.4,
                       'bw': 100,
                       'exp': ['n0s0']})
    b = pds.DataFrame({'tau': 0,
                       'rho': 1,
                       'SNR': -11.8,
                       'bw': 900,
                       'exp': ['n0s0']})
    data = data.append([a, b])
    data['tau'] *= 1e-6
    data['model'] = False

    bw_values = data.bw.unique()

    goal_dprime = 0.770              # equivalent to 70.7% correct

    if not opt:
        dtau = 1 / 24000
        tau_values = np.arange(0, 3e-3 + dtau, dtau)
        # tau_values = np.round(np.arange(0, 3250, 250) * 1e-6, 8)
    else:
        tau_values = data.tau.unique()

    for i_bw, bw in enumerate(bw_values):
        for i_rho, rho in enumerate(rho_values):
            for i_tau, tau in enumerate(tau_values):
                t_param = dict(noise_bw=bw,
                               n_corr=rho,
                               s_ipd=np.pi,
                               itd=tau)
                r_param = dict(noise_bw=bw,
                               n_corr=rho,
                               s_ipd=np.pi,
                               itd=tau)
                snr = find_snr(bin_noise=bin_noise,
                               rho_hat=rho_hat,
                               mon_noise=mon_noise,
                               goal_dprime=goal_dprime,
                               target_param=t_param,
                               reference_param=r_param,
                               ftol=1e-6
                               )
                res_dict = {'tau': tau, 'rho': rho,
                            'SNR': 10 * np.log10(snr), 'bw': bw,
                            'model': True, 'exp': ['n0srho'],
                            'rho_hat': rho_hat, 'mon_noise': mon_noise,
                            'bin_noise': bin_noise}
                data = data.append(pds.DataFrame(res_dict))
        t_param = dict(noise_bw=bw,
                       n_corr=1,
                       s_ipd=0,
                       itd=0)
        r_param = dict(noise_bw=bw,
                       n_corr=1,
                       s_ipd=np.pi,
                       itd=0)
        snr = find_snr(bin_noise=bin_noise,
                       rho_hat=rho_hat,
                       mon_noise=mon_noise,
                       goal_dprime=goal_dprime,
                       target_param=t_param,
                       reference_param=r_param,
                       ftol=1e-6
                       )
        res_dict = {'tau': 0, 'rho': 1,
                    'SNR': 10 * np.log10(snr), 'bw': bw,
                    'model': True, 'exp': ['n0s0'],
                    'rho_hat': rho_hat, 'mon_noise': mon_noise,
                    'bin_noise': bin_noise}
        data = data.append(pds.DataFrame(res_dict))
    if save:
        fname = 'data.h5' if save is True else save
        data.to_hdf(fname, 'bernstein2020')
    return data
