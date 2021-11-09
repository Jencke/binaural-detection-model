"""Calculate model for Rabiner1966."""
from ..model import find_snr
import numpy as np
import pandas as pds
import os


def calc(rho_hat=0.85, bin_noise=0.24, mon_noise=0.71, calc=True,
         save=True, opt=False):
    """Run the experiment."""
    if not calc:
        fname = 'data.h5' if save is True else save
        all_data = pds.read_hdf(fname, 'rabiner1966')
        return all_data

    path, _ = os.path.split(__file__)
    file = os.path.join(path, 'exp_data.csv')
    exp_data = pds.read_csv(file, header=0, skipinitialspace=True)
    exp_data['tau'] = np.round(exp_data['tau'] * 1e-3, 5)
    exp_data['model'] = False
    all_data = exp_data.copy()


    bw = 1100
    fc = 500

    # data = load('fig5', 'b').reset_index()
    # all_data = pds.DataFrame()
    # for i, df in data.iterrows():
    #     res_dict = {'tau': np.round(df.x * 1e-3, 5),
    #                 'BMLD': df.Curve1,
    #                 'model': False}
    #     all_data = all_data.append(pds.DataFrame([res_dict]))

    goal_dprime = 1.5
    #
    dtau = 0.1
    if not opt:
        tau_vals = np.round(np.arange(0, 7.6 + dtau, dtau) * 1e-3, 5)
    else:
        tau_vals = all_data.tau.unique()
    for i_tau, tau in enumerate(tau_vals):
        t_param = dict(noise_bw=bw, n_itd=tau, n_ipd=tau * fc * 2 * np.pi,
                       s_ipd=np.pi)
        r_param = dict(noise_bw=bw, n_itd=tau, n_ipd=tau * fc * 2 * np.pi,
                       s_ipd=np.pi)
        snr = find_snr(bin_noise=bin_noise,
                       rho_hat=rho_hat,
                       mon_noise=mon_noise,
                       goal_dprime=goal_dprime,
                       target_param=t_param,
                       reference_param=r_param)
        res_dict = dict(tau=tau, SNR=10*np.log10(snr[0]), model=True,
                        bin_noise=bin_noise, mon_noise=mon_noise,
                        rho_hat=rho_hat)
        all_data = all_data.append(pds.DataFrame([res_dict]))
    t_param = dict(noise_bw=bw, n_itd=0, s_ipd=0)
    r_param = dict(noise_bw=bw, n_itd=0, s_ipd=0)
    snr_n0s0 = find_snr(bin_noise=bin_noise,
                        rho_hat=rho_hat,
                        mon_noise=mon_noise,
                        goal_dprime=goal_dprime,
                        target_param=t_param,
                        reference_param=r_param)

    exp_data = all_data[~all_data.model]
    snr = 10*np.log10(snr_n0s0[0]) - exp_data.BMLD
    all_data.loc[~all_data.model, 'SNR'] = snr
    all_data = all_data.drop(columns='BMLD')
    if save:
        fname = 'data.h5' if save is True else save
        all_data.to_hdf(fname, 'rabiner1966')
    return all_data
