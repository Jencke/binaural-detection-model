import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
import numpy as np
from  .calc import calc
from ..model_helpers import create_figure

def plot(gridspec=None, file='data.h5', marker='o'):
    if gridspec is None:
        gridspec = create_figure()

    all_data = calc(calc=False, save=file)

    model_data = all_data[all_data.model]
    paper_data = all_data[~all_data.model]

    colors = ['#727272',
              '#599ad3',
              '#79c36a',
              '#f9a65a',
              '#f1595f',
              '#9e66ab',
              '#cd7058',
              '#d77fb3']

    ax = plt.subplot(gridspec)
    for i, (bw, g) in enumerate(model_data.groupby('bw')):
        ax.plot(g.rho, g.SNR, color=colors[i],
                label=f"{bw} Hz")
    for i, (bw, g) in enumerate(paper_data.groupby('bw')):
        ax.plot(g.rho, g.SNR, marker=marker, color=colors[i],
                ls='')
    ax.set_ylim(-30, 0)
    ax.set_xlim(-1.05, 1.05)
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax.yaxis.set_major_locator(plt.MultipleLocator(5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(2.5))
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel('SNR / dB')
    ax.legend(ncol=4, handlelength=0.8, columnspacing=0.8,
              loc =(-0.27, 1), handletextpad=0.5)
    ax.set_title('Bernstein \& Trahiotis, 2014', pad=15)
    return ax
