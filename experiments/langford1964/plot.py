import matplotlib.pyplot as plt
from .calc import calc
from ..model_helpers import create_figure

def plot(gridspec=None, file='data.h5', marker='o'):
    if gridspec is None:
        gridspec = create_figure()
    all_data = calc(calc=False, save=file)
    ax = plt.subplot(gridspec)
    for i_exp, (exp, exp_g) in enumerate(all_data.groupby('exp')):
        mod_data = exp_g[exp_g.model]
        exp_data = exp_g[~exp_g.model]
        p  = ax.plot(mod_data.tau * 1e3, mod_data.SNR, '-', label=exp)
        ax.plot(exp_data.tau*1e3, exp_data.SNR, marker=marker,
                color=p[0].get_color(), ls='')
    ax.set_xlim(-0.1, 9.1)
    ax.set_ylim(-25, -9.5)
    ax.yaxis.set_major_locator(plt.MultipleLocator(5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(2.5))
    ax.xaxis.set_major_locator(plt.MultipleLocator(2))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.set_xlabel(r'$\Delta t$ / ms')
    ax.set_ylabel('SNR / dB')
    ax.legend(ax.get_lines()[::2], [r'$N_{\Delta t} S_0$', r'$N_{\Delta t} S_\pi$'],
              handlelength=1, loc=(0.15, 1), ncol=3,
              columnspacing=1)
    ax.set_title('Langford \& Jeffress, 1964', pad=15)
    return ax
