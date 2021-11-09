"""Calculate model for van der Heijden1999."""
from ..model import find_snr
import numpy as np
import os
import pandas as pds
import audiotools as audio


def calc(rho_hat=0.90, bin_noise=0.19, mon_noise=0.61, calc=True,
                   save=True):
    """Run the experiment."""
    if not calc:
        fname = 'data.h5' if save is True else save
        all_data = pds.read_hdf(fname, 'vanderhejden1999')
        return all_data

    path, _ = os.path.split(__file__)
    file = os.path.join(path, 'exp_data.csv')
    exp_data = pds.read_csv(file, header=0, skipinitialspace=True)

    bw = 900
    n_band_level = 45.5
    goal_dprime = 0.770              # equivalent to 70.7% correct
    # tau_max = 4.2e-3
    # tau_vals = np.linspace(0, tau_max, 51)
    tau_vals = exp_data['tau'].values * 1e-6
    ftol = 1e-3

    snr_n_tau_s_pi = np.zeros(len(tau_vals))
    snr_n_tau_s_0 = np.zeros(len(tau_vals))
    ddn_s_pi = np.zeros(len(tau_vals))
    ddn_s_0 = np.zeros(len(tau_vals))
    snr_n_0_s_tau = np.zeros(len(tau_vals))
    snr_n_pi_s_tau = np.zeros(len(tau_vals))
    for i_tau, tau in enumerate(tau_vals):
        t_param = dict(noise_bw=bw, n_itd=tau, s_ipd=np.pi)
        r_param = dict(noise_bw=bw, n_itd=tau, s_ipd=np.pi)
        snr_n_tau_s_pi[i_tau] = find_snr(bin_noise=bin_noise,
                                         rho_hat=rho_hat,
                                         mon_noise=mon_noise,
                                         goal_dprime=goal_dprime,
                                         target_param=t_param,
                                         reference_param=r_param,
                                         ftol=ftol)
        t_param = dict(noise_bw=bw, n_itd=tau, s_ipd=0)
        r_param = dict(noise_bw=bw, n_itd=tau, s_ipd=0)
        snr_n_tau_s_0[i_tau] = find_snr(bin_noise=bin_noise,
                                        rho_hat=rho_hat,
                                        mon_noise=mon_noise,
                                        goal_dprime=goal_dprime,
                                        target_param=t_param,
                                        reference_param=r_param,
                                        ftol=ftol)
        t_param = dict(noise_bw=bw, n_itd=tau, s_ipd=np.pi, ddn=True)
        r_param = dict(noise_bw=bw, n_itd=tau, s_ipd=np.pi, ddn=True)
        ddn_s_pi[i_tau] = find_snr(bin_noise=bin_noise,
                                             rho_hat=rho_hat,
                                             mon_noise=mon_noise,
                                             goal_dprime=goal_dprime,
                                             target_param=t_param,
                                             reference_param=r_param,
                                             ftol=ftol)
        t_param = dict(noise_bw=bw, n_itd=tau, s_ipd=0, ddn=True)
        r_param = dict(noise_bw=bw, n_itd=tau, s_ipd=0, ddn=True)
        ddn_s_0[i_tau] = find_snr(bin_noise=bin_noise,
                                            rho_hat=rho_hat,
                                            mon_noise=mon_noise,
                                            goal_dprime=goal_dprime,
                                            target_param=t_param,
                                            reference_param=r_param,
                                            ftol=ftol)
        t_param = dict(noise_bw=bw, s_itd=tau, s_ipd=np.pi)
        r_param = dict(noise_bw=bw, s_itd=tau, s_ipd=np.pi)
        snr_n_0_s_tau[i_tau] = find_snr(bin_noise=bin_noise,
                                        rho_hat=rho_hat,
                                        mon_noise=mon_noise,
                                        goal_dprime=goal_dprime,
                                        target_param=t_param,
                                        reference_param=r_param,
                                        ftol=ftol)
        t_param = dict(noise_bw=bw, s_itd=tau, n_ipd=np.pi)
        r_param = dict(noise_bw=bw, s_itd=tau, n_ipd=np.pi)
        snr_n_pi_s_tau[i_tau] = find_snr(bin_noise=bin_noise,
                                         rho_hat=rho_hat,
                                         mon_noise=mon_noise,
                                         goal_dprime=goal_dprime,
                                         target_param=t_param,
                                         reference_param=r_param,
                                         ftol=ftol)
    noise_level = audio.band2rms(n_band_level, bw)

    n_tau_s_0_model = pds.DataFrame({'tau': tau_vals,
                                     'SNR': 10 * np.log10(snr_n_tau_s_0),
                                     'exp': 'n_tau_s_0',
                                     'model': True,
                                     'rho_hat': rho_hat,
                                     'mon_noise': mon_noise,
                                     'bin_noise': bin_noise
                                     })
    ddn_s_0_model = pds.DataFrame({'tau': tau_vals,
                                   'SNR': 10 * np.log10(ddn_s_0),
                                   'exp': 'ddn_s_0',
                                   'model': True,
                                   'rho_hat': rho_hat,
                                   'mon_noise': mon_noise,
                                   'bin_noise': bin_noise
                                   })
    n_pi_s_tau_model = pds.DataFrame({'tau': tau_vals,
                                      'SNR': 10 * np.log10(snr_n_pi_s_tau),
                                      'exp': 'n_pi_s_tau',
                                      'model': True,
                                      'rho_hat': rho_hat,
                                      'mon_noise': mon_noise,
                                      'bin_noise': bin_noise
                                      })
    n_pi_s_tau_data = pds.DataFrame({'tau': exp_data.tau * 1e-6,
                                     'SNR': exp_data['NpiStau'] - noise_level,
                                     'exp': 'n_pi_s_tau',
                                     'model': False
                                     })
    n_tau_s_0_data = pds.DataFrame({'tau': exp_data.tau * 1e-6,
                                    'SNR': exp_data['NtauS0'] - noise_level,
                                    'exp': 'n_tau_s_0',
                                    'model': False
                                    })
    n_tau_s_pi_model = pds.DataFrame({'tau': tau_vals,
                                      'SNR': 10 * np.log10(snr_n_tau_s_pi),
                                      'exp': 'n_tau_s_pi',
                                      'model': True,
                                      'rho_hat': rho_hat,
                                      'mon_noise': mon_noise,
                                      'bin_noise': bin_noise
                                      })
    ddn_s_pi_model = pds.DataFrame({'tau': tau_vals,
                                    'SNR': 10 * np.log10(ddn_s_pi),
                                    'exp': 'ddn_s_pi',
                                    'model': True,
                                    'rho_hat': rho_hat,
                                    'mon_noise': mon_noise,
                                    'bin_noise': bin_noise
                                    })
    n_tau_s_pi_data = pds.DataFrame({'tau': exp_data.tau * 1e-6,
                                     'SNR': exp_data['NtauSpi'] - noise_level,
                                     'exp': 'n_tau_s_pi',
                                     'model': False
                                     })
    n_0_s_tau_model = pds.DataFrame({'tau': tau_vals,
                                     'SNR': 10 * np.log10(snr_n_0_s_tau),
                                     'exp': 'n_0_s_tau',
                                     'model': True,
                                     'rho_hat': rho_hat,
                                     'mon_noise': mon_noise,
                                     'bin_noise': bin_noise
                                     })
    n_0_s_tau_data = pds.DataFrame({'tau': exp_data.tau * 1e-6,
                                    'SNR': exp_data['N0Stau'] - noise_level,
                                    'exp': 'n_0_s_tau',
                                    'model': False
                                    })
    all_data = pds.DataFrame().append([n_tau_s_0_model,
                                       n_tau_s_pi_model,
                                       n_tau_s_0_data,
                                       n_tau_s_pi_data,
                                       ddn_s_pi_model,
                                       ddn_s_0_model,
                                       n_0_s_tau_model,
                                       n_pi_s_tau_model,
                                       n_pi_s_tau_data,                                                                          n_0_s_tau_data
                                       ])
    all_data['bw'] = len(all_data) * [bw]
    # all_data['paper'] = len(all_data) * ['vanderhejden1999']

    if save:
        fname = 'data.h5' if save is True else save
        all_data.to_hdf(fname, 'vanderhejden1999')

    return all_data
