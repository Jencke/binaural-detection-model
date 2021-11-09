"""Complex cross-correlation based BMLD Model."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from .model_helpers import calc_bandwidth


def phase_spectrum(freq, itd=0, ipd=0):
    """Phase sepectrum for a given ITD and IPD.

    Paramters
    ---------
    freq : ndarray
        The frequency in Hz
    itd : scalar
        The time difference between the two channels
    ipd : scalar
        The frequency independent phase shift

    Returns
    -------
    The phase spectrum in rad : ndarray
    """
    omega = 2 * np.pi * freq
    ipd = ipd + omega * itd
    # wrap the ipd
    ipd = (ipd + np.pi) % (2 * np.pi) - np.pi
    return ipd


def cmplx_corr_coeff(noise_bw, snr=0, itd=0, ipd=0, s_ipd=np.pi, s_itd=0,
                     n_itd=0, n_ipd=0, n_corr=1, fc=500, ddn=False, **kwargs):
    """Complex correlation coefficent for tone in noise experiments.

    Parameters
    ----------
    noise_bw : scalar
        bandwidth of the rectangular noise in Hz
    snr : scalar
        Signal to noise Ratio (not in dB)
    itd : scalar
        ITD in seconds for both noise and tone (default=0)
    ipd : scalar
        IPD in rad for both noise and tone (default=0)
    s_itd : scalar
        ITD for the tone in seconds. (default=0)
    s_ipd : scalar
        IPD for the tone in rad. (default=pi)
    n_itd : scalar
        ITD for the noise in seconds. (default=0)
    n_ipd : scalar
        IPD for the noise in seconds. (default=0)
    n_corr : scalar
        Correlation of the noise. (default = 1)
    fc : scalar
        center frequency of the noise and peripheral filter
    ddn: bool
        Indicates if to calculate results for double delayed noise,
        (default=False)

    Returns
    -------
    The complex correlation-coefficent : complex

    """
    f_low = fc - noise_bw / 2
    f_high = fc + noise_bw / 2

    # The tone PSD has thethe SNR with the IPD
    tone_itd = -(s_itd + itd)
    tone_ipd = s_ipd + ipd
    tone_phase = phase_spectrum(fc, itd=tone_itd, ipd=tone_ipd)
    tone_psd = snr * np.exp(1j * tone_phase)

    # Integrate over the noise psd. this is done separately for real and
    # imaginary part
    noise_itd = -(itd + n_itd)     # overall noise itd
    noise_ipd = n_ipd + ipd        # overall noise ipd
    noise_psd_func = _noise_psd(fc, itd=noise_itd, ipd=noise_ipd, ddn=ddn)

    def real_func(freq):
        return np.real(noise_psd_func(freq))

    def imag_func(freq):
        return np.imag(noise_psd_func(freq))

    real_int = quad(real_func, f_low, f_high)[0]
    imag_int = quad(imag_func, f_low, f_high)[0]
    noise_int = real_int + 1j * imag_int
    # the noise correlation has to be normalized by the bandwidth to make sure
    # that the overall noise energy equaled 1. It is also multiplied by the
    # noise correlation as only coherent noise contibutes to the CSPD
    noise_corr = n_corr * noise_int / noise_bw

    # add the tone_psd to the noise integral to get the non-normalized CCC
    ccc_num = noise_corr + tone_psd

    # calculate the noise power
    mon_psd = _noise_psd(fc, itd=0, ipd=0)

    def mon_func(freq):
        return np.real(mon_psd(freq))

    p_noise = quad(mon_func, f_low, f_high)[0] / noise_bw

    # Normalize the CCC numerator by the total power
    p_total = p_noise + snr
    ccc = ccc_num / (p_total)

    return ccc


def calc_gtf_gain(freq, f_0, bw, n=4):
    """Calculate the gain of a gammatone filter.

    Parameters
    ----------
    freq : float
      Frequency in Hz
    f_0 : float
      the center frequency
    bw : float
      Equivalent rectangular bandwidth in Hz
    n : int
      Filter order (default = 4)
    """
    alpha = ((np.pi * np.math.factorial(2 * n - 2) * 2**(2 - 2 * n))
             / np.math.factorial(n - 1)**2)

    b = bw / alpha  # calculate b for bw=equivalent-rectangular bandwidth

    gain = (1 + (freq - f_0)**2 / b**2)**(-n/2)
    return gain


def calc_eff_snr(noise_bw, snr, fc=500):
    """Effective power of a noiseband after peripheral filtering.

    Parameters
    ----------
    noise_bw: scalar
        bandwidth of the noise
    fc : scalar
        Center frequency of the noiseband
    bandlevel : scalar
        Bandlevel of the noise

    Returns
    -------
    scalar : The effective stimulus power
    """
    # calculate the amplitude of the noiseband
    # noise energy would be amplitude squared times bandwidth

    gtf_bw = calc_bandwidth(fc, 'erb')

    # Due to the gammatone filter, the "effective" bandwidth is
    # reduced. in order of calculating the effective bandwidth, it is
    # necessary to integrate over the gammatone filters pdf. For large
    # noise bandwidth, the effective bandwidth converges towards the
    # ERB bandwidth of the filter
    def int_func(freq):
        return calc_gtf_gain(freq, fc, gtf_bw)**2

    eff_bw = quad(int_func, fc - noise_bw/2, fc+noise_bw/2)[0]

    # Power is calculated:
    eff_snr = snr * noise_bw / eff_bw
    return eff_snr


def calc_mon_dprime(mon_noise, param):
    """Calculate the monaural dprime"""
    snr = param['snr']
    noise_bw = param['noise_bw']
    fc = 500 if 'fc' not in param else param['fc']

    eff_snr = calc_eff_snr(noise_bw, snr=snr, fc=fc)
    d_mon = eff_snr / mon_noise
    return d_mon


def calc_bin_dprime(rho_hat, bin_noise, target_param, reference_param):
    """Calculate the binaural dprime."""
    # calculate CCC for reference and target
    t_ccc = cmplx_corr_coeff(**target_param)
    r_ccc = cmplx_corr_coeff(**reference_param)
    # calculate absolute value and argument of the CCC
    t_angle = np.angle(t_ccc)
    r_angle = np.angle(r_ccc)
    t_abs = np.abs(t_ccc)
    r_abs = np.abs(r_ccc)
    # apply z-transform to the absolute_values and multiply with the
    # binaural noise which limits maximal performance, then add the
    # orignal phase angle
    t_ztrafo = np.arctanh(t_abs * rho_hat) * np.exp(1j * t_angle)
    r_ztrafo = np.arctanh(r_abs * rho_hat) * np.exp(1j * r_angle)

    # dprime is given by the eucledian distance in the complex plane
    # divided by the deciscion noise
    d_bin = np.abs(t_ztrafo - r_ztrafo) / bin_noise

    return d_bin


def calc_dprime(bin_noise, rho_hat, mon_noise,
                target_param,
                reference_param):
    """Calculate overall sensitivity index."""
    d_bin = calc_bin_dprime(rho_hat=rho_hat,
                            bin_noise=bin_noise,
                            target_param=target_param,
                            reference_param=reference_param)
    d_mon = calc_mon_dprime(mon_noise, target_param)
    d_sum = np.sqrt(d_mon**2 + d_bin**2)
    return d_sum


def find_snr(bin_noise, rho_hat, mon_noise, goal_dprime,
             target_param, reference_param, ftol=1e-6):
    """Find the SNR at which the overall d' reached the goal_dprime.

    Parameters
    ----------
    bin_noise : scalar
        Standard deviation of the noise that determines binaural sensitivity
    rho_hat : scalar
        The factor in the range [0, 1] that determines maximal sensitivy to
        changes in coherence.
    mon_noise : scalar
        Standard deviation of the noise that determines monaural sensitivity
    goal_dprime : scalar
        The sensitivity index d' for which to determine the snr
    target_param : dict
        A dictionary containing parameters for the cmplx_corr_coeff function
        which define the target stimulus
    reference_param : dict
        A dictionary containing parameters for the cmplx_corr_coeff function
        which define the target stimulus
    ftol : scalar, optional
        The tolarance of the numeric optimizer (default=1e-6)

    Returns:
    --------
    The determined SNR : scalar

    """
    def to_minimize(snr):
        target_param['snr'] = 10**(snr / 10)
        reference_param['snr'] = 0
        dprime = calc_dprime(bin_noise=bin_noise,
                             rho_hat=rho_hat,
                             mon_noise=mon_noise,
                             target_param=target_param,
                             reference_param=reference_param)
        diff = np.abs(dprime - goal_dprime)
        # print()
        # print(f'{diff=}')
        # print(f'{snr=}')
        return diff

    p = [10*np.log10(0.01)]
    bounds = [[10*np.log10(1e-6), 10*np.log10(mon_noise)]]
    snr = minimize(to_minimize,
                   p,
                   bounds=bounds,
                   method='SLSQP',
                   options=dict(ftol=ftol))

    if not snr.success:
        print('optim. problem')

    return 10**(snr.x/10)


def find_rho_jnd(bin_noise, rho_hat, mon_noise, goal_dprime,
             target_param, reference_param, ftol=1e-6):
    """Find the correlation detection JND.

    Parameters
    ----------
    bin_noise : scalar
        Standard deviation of the noise that determines binaural sensitivity
    rho_hat : scalar
        The factor in the range [0, 1] that determines maximal sensitivy to
        changes in coherence.
    mon_noise : scalar
        Standard deviation of the noise that determines monaural sensitivity
    goal_dprime : scalar
        The sensitivity index d' for which to determine the snr
    target_param : dict
        A dictionary containing parameters for the cmplx_corr_coeff function
        which define the target stimulus
    reference_param : dict
        A dictionary containing parameters for the cmplx_corr_coeff function
        which define the target stimulus
    ftol : scalar, optional
        The tolarance of the numeric optimizer (default=1e-6)

    Returns:
    --------
    The determined JND : scalar

    """
    def to_minimize(delta_corr):
        target_param['n_corr'] = reference_param['n_corr'] + delta_corr
        dprime = calc_dprime(bin_noise=bin_noise,
                             rho_hat=rho_hat,
                             mon_noise=mon_noise,
                             target_param=target_param,
                             reference_param=reference_param)
        diff = np.abs(dprime - goal_dprime)
        return diff

    p = [0.5 * (1 - target_param['n_corr'])]
    bounds = [[0, 1 - target_param['n_corr']]]
    jnd = minimize(to_minimize,
                   p,
                   bounds=bounds,
                   method='SLSQP',
                   options=dict(ftol=ftol))
    # print(jnd)
    if not jnd.success:
        print('optim. problem')

    return jnd.x


def _noise_psd(fc, itd, ipd, ddn=False):
    gtf_bw = calc_bandwidth(fc, 'erb')
    # gammatone_filter as modiolus, ipd as phase
    if not ddn:
        def func(freq):
            return (calc_gtf_gain(freq, fc, gtf_bw)**2 *
                    np.exp(1j * phase_spectrum(freq, itd, ipd)))
    else:
        def func(freq):
            return (0.5 * calc_gtf_gain(freq, fc, gtf_bw)**2 *
                    np.exp(1j * phase_spectrum(freq, itd, ipd))
                    +
                    0.5 * calc_gtf_gain(freq, fc, gtf_bw)**2 *
                    np.exp(1j * phase_spectrum(freq, -itd, ipd)))
    return func
