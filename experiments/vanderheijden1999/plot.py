import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
from .calc import calc
from ..model_helpers import create_figure


def plot(gridspec=None, file='data.h5', marker='o'):
    if gridspec is None:
        gridspec = create_figure()

    all_data = calc(calc=False, save=file)

    mod_data = all_data[all_data.model]
    exp_data = all_data[~all_data.model]
    #
    mod_tau = all_data[all_data.model].tau.unique() * 1e3
    exp_tau = all_data[~all_data.model].tau.unique() * 1e3


    gs = GridSpecFromSubplotSpec(2, 1, gridspec)
    ax = [plt.subplot(gs[0]), plt.subplot(gs[1])]
    # fig, ax = plt.subplots(2, 1, sharex='all', sharey='all')
    p = ax[0].plot(mod_tau, mod_data[mod_data.exp == 'n_tau_s_pi'].SNR,
                   label=r'$N_{\Delta t} S_\pi$')
    ax[0].plot(exp_tau, exp_data[exp_data.exp == 'n_tau_s_pi'].SNR,
               marker=marker, color=p[0].get_color())
    p = ax[0].plot(mod_tau, mod_data[mod_data.exp == 'n_0_s_tau'].SNR,
                   label=r'$N_{\pi} S_{\Delta t}$')
    ax[0].plot(exp_tau, exp_data[exp_data.exp == 'n_0_s_tau'].SNR,
               marker=marker, color=p[0].get_color())
    p = ax[1].plot(mod_tau, mod_data[mod_data.exp == 'n_tau_s_0'].SNR,
                   label=r'$N_{\Delta t} S_0$')
    ax[1].plot(exp_tau, exp_data[exp_data.exp == 'n_tau_s_0'].SNR,
               marker=marker, color=p[0].get_color())
    p = ax[1].plot(mod_tau, mod_data[mod_data.exp == 'n_pi_s_tau'].SNR,
                   label=r'$N_{0} S_{\Delta t}$')
    ax[1].plot(exp_tau, exp_data[exp_data.exp == 'n_pi_s_tau'].SNR,
               marker=marker, color=p[0].get_color())
    for a in ax:
        a.set_xlim(-0.05, 4.25)
        a.set_ylim(-30, -12.5)
        a.yaxis.set_major_locator(plt.MultipleLocator(5))
        a.yaxis.set_minor_locator(plt.MultipleLocator(2.5))
        a.xaxis.set_major_locator(plt.MultipleLocator(1))
        a.set_ylabel('SNR / dB')
    ax[0].set_xticklabels([])
    ax[1].set_xlabel('$\Delta t$ / ms')
    ax[0].set_title('van der Heijden \& Trahiotis 1999', pad=15)
    # ax[0].legend(handlelength=1, loc=(0.07, 0.95), ncol=2)
    # ax[1].legend(handlelength=1, loc=(0.07, 0.9), ncol=2)
    ax[0].legend()
    ax[1].legend()
    return ax
