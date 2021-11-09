import matplotlib.pyplot as plt
from .calc import calc
from ..model_helpers import create_figure


def plot(gridspec=None, file='data.h5', marker='o'):
    if gridspec is None:
        gridspec = create_figure()

    all_data = calc(calc=False, save=file)
    ax = plt.subplot(gridspec)
    for dprime, dprime_g in all_data.groupby('dprime'):
        mod_data = dprime_g[dprime_g.model]
        exp_data = dprime_g[~dprime_g.model]
        jnd = (mod_data['corr'] + mod_data['jnd'])**2 - mod_data['corr']**2
        p=plt.plot(mod_data['corr'], mod_data.jnd, label=rf"$d'={dprime}$")
        # p = plt.plot(mod_data['corr'], jnd)
        plt.plot(exp_data['corr'], exp_data.jnd, marker=marker,
                 ls='', color=p[0].get_color())
    # ax.legend(handlelength=1, ncol=3, )
    ax.legend(handlelength=1, loc=(-0.25, 1), ncol=3,
              columnspacing=1)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(0.02, 1)
    ax.set_yscale('log')
    ax.set_xlabel(r'$\rho_{ref}$')
    ax.set_ylabel(r'$\Delta\rho$')
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # ax.grid(which='both', ls='-', color='#DDDDDD')
    ax.set_yticks([0.03, 0.1, 0.30, 1])
    ax.set_yticklabels([0.03, 0.1, 0.30, 1])
    ax.set_title(r'Pollack \& Trittipoe 1959', pad=15)
    return ax
