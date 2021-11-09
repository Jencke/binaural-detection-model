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
    for i_bw, (bw, bw_g) in enumerate(all_data.groupby('bw')):
        for i_rho, (rho, rho_g) in enumerate(bw_g[bw_g.exp=='n0srho'].groupby('rho')):
            mod_data = rho_g[rho_g.model]
            exp_data = rho_g[~rho_g.model]
            p = ax[i_bw].plot(mod_data.tau * 1e3, mod_data.SNR,
                              label=rf'$\rho$={rho}')
            ax[i_bw].plot(exp_data.tau * 1e3, exp_data.SNR,
                          marker=marker, ls='',
                          color=p[0].get_color())
        n0s0=bw_g[bw_g.exp=='n0s0']
        ax[i_bw].hlines([n0s0[n0s0.model].SNR], 0, 3, color='k')
        ax[i_bw].plot(exp_data.tau * 1e3, len(exp_data.tau) * [n0s0[~n0s0.model].SNR],
                      marker=marker, color='k', ls='')
    for a in ax:
        a.set_xlim(-0.05, 3.05)
        a.yaxis.set_major_locator(plt.MultipleLocator(5))
        a.yaxis.set_minor_locator(plt.MultipleLocator(2.5))
        a.xaxis.set_major_locator(plt.MultipleLocator(1))
        a.set_ylabel('SNR / dB')
    ax[0].set_ylim(-20, 0)
    ax[1].set_ylim(-28, -8)
    for i_bw, bw in enumerate(all_data['bw'].unique()):
        a = ax[i_bw]
        y_lim = a.get_ylim()
        y_range = y_lim[1] - y_lim[0]
        t_pos = y_lim[0] + y_range * 0.05
        a.text(2.5, t_pos, f'{bw}\,Hz', fontsize=8)
    ax[0].set_xticklabels([])
    ax[1].set_xlabel('$\Delta t$ / ms')
    ax[0].set_title('Bernstein \& Trahiotis 2020', pad=15)
    ax[0].legend(ncol=4, handlelength=0.8, columnspacing=0.8,
              loc =(-0.27, 1), handletextpad=0.5)
    # ax[1].legend(handlelength=1, loc=(0.07, 0.9), ncol=2)
    return ax
