import matplotlib.pyplot as plt
from .calc import calc
from ..model_helpers import create_figure


def plot(gridspec=None, file='data.h5', marker='o'):
    if gridspec is None:
        gridspec = create_figure()
    all_data = calc(calc=False, save=file)
    mod_data = all_data[all_data.model]
    exp_data = all_data[~all_data.model]
    ax = plt.subplot(gridspec)
    # shift = (exp_data[exp_data.tau==0].BMLD[0]
    #          + mod_data[mod_data.tau==0].SNR[0])
    p = ax.plot(mod_data.tau*1e3, mod_data.SNR)
    ax.plot(exp_data.tau*1e3, exp_data.SNR, marker=marker,
            color=p[0].get_color(), ls='')
    ax.set_xlim(-0.1, 8.0)
    ax.set_ylim(-25, -15)
    ax.yaxis.set_major_locator(plt.MultipleLocator(5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(2.5))
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.set_ylabel('SNR / dB')
    ax.set_xlabel(r'$\tau$ / ms')
    ax.set_title('Rabiner et al., 1966 ', pad=15)
    return ax
