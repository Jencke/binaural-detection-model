"""Some helper functions for the model."""
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec

def calc_bandwidth(fc, scale='cbw'):
    r"""Calculate approximation of auditory filter bandwidth.

    This Function calculates aproximations for the auditory filter
    bandwidth using differnt concepts:

     - cbw: Use the critical bandwidth concept following [1]_
     - erb: Use the equivalent rectangular bandwith concept following
       [2]_

    Equation used for critical bandwidth:
    .. math:: B = 25 + 75 (1 + 1.4 \frac{f_c}{1000}^2)^0.69

    Equation used for critical equivalent rectangular bandwith:
    .. math:: B = 24.7 (4.37 \frac{f_c}{1000} + 1)

    Parameters
    -----------
    fc : float or ndarray
      center frequency in Hz

    scale : str
      String indicating the scale that should be used possible values:
      'cbw' or 'erb'. (default='cbw')

        ..[1] Zwicker, E., & Terhardt, E. (1980). Analytical
              expressions for critical-band rate and critical
              bandwidth as a function of frequency. The Journal of the
              Acoustical Society of America, 68(5), 1523-1525.

        ..[2] Glasberg, B. R., & Moore, B. C. (1990). Derivation of
              auditory filter shapes from notched-noise data. Hearing
              Research, 47(1-2), 103-138.

    """
    if 'cbw' in scale:
        bw = 25 + 75 * (1 + 1.4 * (fc / 1000.)**2)**0.69
    elif 'erb' in scale:
        bw = 24.7 * (4.37 * (fc / 1000.) + 1)

    return bw

def create_figure():
    fig = plt.figure(figsize=(4, 3))
    main_gs = plt.GridSpec(1, 1, left=0.2, right=0.98, bottom=0.22, top=0.82)
    gs = GridSpecFromSubplotSpec(1, 1, main_gs[0])
    return gs[0]
