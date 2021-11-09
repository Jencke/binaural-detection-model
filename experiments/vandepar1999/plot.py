import matplotlib.pyplot as plt
from .calc import calc
import numpy as np
from ..model_helpers import create_figure


def plot(gridspec=None, file='data.h5', marker='o'):
    if gridspec is None:
        gridspec = create_figure()
    all_data = calc(calc=False, save=file)
    ax = plt.subplot(gridspec)
    for exp, exp_g in all_data.groupby('exp'):
        mod_data = exp_g[exp_g.model]
        exp_data = exp_g[~exp_g.model]
        p = ax.plot(mod_data.bw, mod_data.SNR)
        ax.plot(exp_data.bw, exp_data.SNR, marker=marker,
                ls='', color=p[0].get_color())
    f = np.logspace(np.log10(5), np.log10(100))
    ax.plot(f, -1.5 * np.log2(f) + 6.5, ':', color='k')
    f = np.logspace(np.log10(100), np.log10(1000))
    ax.plot(f, -3 * np.log2(f) + 16.5, ':', color='k')
    ax.text(15, -1, r'-1.5\,dB/oct', fontsize=6, rotation=-11)
    ax.text(200, -9.4, r'-3\,dB/oct', fontsize=6, rotation=-20)
    ax.legend(np.asarray(ax.lines)[[2, 4, 0]],
              [r'$N_0 S_{\pi}$', r'$N_{\pi} S_{0}$', r'$N_{0} S_{0}$'],
              handlelength=1, ncol=3, loc=(-0.08, 1), columnspacing=1)
    ax.set_xscale('log')
    ax.set_ylim(-30, 5)
    ax.yaxis.set_major_locator(plt.MultipleLocator(5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(2.5))
    ax.xaxis.set_major_locator(plt.LogLocator(base=10, subs=(1, ), numticks=10))
    ax.set_xlim(4, 1100)
    ax.set_title('van de Par \& Kohlrausch 1999', pad=15)
    ax.set_xlabel('Bandwidth / Hz')
    ax.set_ylabel('SNR / dB')
    return ax

if __name__ == '__main__':

    fig = plt.figure(figsize=(2, 1.5))
    main_gs = plt.GridSpec(1, 1, left=0.2, right=0.98, bottom=0.22, top=0.82)
    gs = GridSpecFromSubplotSpec(1, 1, main_gs[0])
    plot_results(gs[0])
