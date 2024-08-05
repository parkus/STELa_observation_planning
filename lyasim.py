from pathlib import Path

import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy import table
from scipy.special import wofz
from math import pi, cos, sin

import utilities as utils

reffolder = Path('reference_tables')

# region intrinsic lya line shape
fwhm_broad = 400*u.km/u.s # rough average from youngblood+ 2016
fwhm_narrow = 125*u.km/u.s
amp_ratio_broad_narrow = 0.04
sig_ratio_reverse_narrow = 0.5 # just a by-eye calibration using the results of Wood+ 2005
amp_ratio_reverse_narrow = 0.3/0.5
#endregion

#region parameters you probably don't care about
wlab_H = 1215.67*u.AA
wlab_D = 1215.34*u.AA
D_H = 1.5e-5  # from bourrier+ 2017 (kepler 444)
f = 4.1641e-01
A = 6.2649e8 / u.s
mH, mD = 1 * u.u, 2 * u.u
#endregion

# region parameters you  might care about
Tism = 10000*u.K
# endregion

#region Lya transit analysis
def ism_velocity(ra, dec):
    # see http://lism.wesleyan.edu/ColoradoLIC.html and
    # lallement & bertin 1992
    if hasattr(ra, '__iter__'):
        raise ValueError('No vectors or columns please.')
    ra *= pi/180
    dec *= pi/180
    v_sun_ra = 74.5*pi/180
    v_sun_dec = -6*pi/180
    sun_v = np.array((cos(v_sun_ra)*cos(v_sun_dec),
                      sin(v_sun_ra)*cos(v_sun_dec),
                      sin(v_sun_dec)))
    star_v = np.array((cos(ra)*cos(dec),
                       sin(ra)*cos(dec),
                       sin(dec)))
    rv = np.dot(sun_v, star_v) * 26.*u.km/u.s
    return rv

def transmission(w, rv, Nh, T):
    w0s = doppler_shift((wlab_H, wlab_D)*u.AA, rv)
    xsections = [voigt_xsection(w, w0, f, A, T, m) for
                 w0, m in zip(w0s, (mH, mD))]
    tau = xsections[0]*Nh + xsections[1]*Nh*D_H
    return np.exp(-tau)


def w2v(w):
    return  (w/wlab_H.value - 1)*const.c.to('km/s').value


def v2w(v):
    return  (v/const.c.to('km/s').value + 1)*wlab_H.value


def reversed_lya_profile(w, rv, flux):
    w0 = doppler_shift(wlab_H, rv)

    def sig_from_fwhm(fwhm):
        sig = fwhm/2/np.sqrt(2*np.log(2))*w0/const.c
        return sig.to(w0.unit)
    sig_broad, sig_narrow = map(sig_from_fwhm, (fwhm_broad, fwhm_narrow))
    sig_reverse = sig_narrow*sig_ratio_reverse_narrow

    def gaussian_profile(amp, sigma):
        f = amp*np.exp(-(w - w0)**2/2/sigma**2)
        return f

    def gaussian_flux(amp, sigma):
        return np.sqrt(2*np.pi)*sigma*amp

    raw_flux = (gaussian_flux(1, sig_narrow)
                + gaussian_flux(amp_ratio_broad_narrow, sig_broad)
                - gaussian_flux(amp_ratio_reverse_narrow, sig_reverse))
    raw_profile = (gaussian_profile(1, sig_narrow)
                   + gaussian_profile(amp_ratio_broad_narrow, sig_broad)
                   - gaussian_profile(amp_ratio_reverse_narrow, sig_reverse))
    profile = flux/raw_flux * raw_profile
    return profile


def lya_at_earth(wgrid, Flya, rv_star, rv_ism, Nh, Tism):
    intrinsic = reversed_lya_profile(wgrid, rv_star, Flya)
    transmitted = transmission(wgrid, rv_ism, Nh, Tism)
    observed = intrinsic * transmitted
    return observed


wgrid_std = np.arange(1210., 1220., 0.005) * u.AA
def lya_at_earth_auto(catalog_row, n_H, lya_factor):
    planet = catalog_row
    Flya = planet['Flya_at_earth'] * u.Unit('erg s-1 cm-2') * lya_factor
    rv_star = planet['st_radv'] * u.km / u.s
    rv_ism = ism_velocity(planet['ra'], planet['dec'])
    Nh = n_H * planet['sy_dist'] * u.pc
    observed = lya_at_earth(wgrid_std, Flya, rv_star, rv_ism, Nh, Tism)
    return observed


def depth_extended_hill_transit(planet_row):
    s = planet_row
    Mp = s['pl_bmasse'] * u.Mearth
    if Mp > 0.41*u.Mjup and s['pl_bmasse'] not in ['Mass', 'Msini']: # based on Chen and Kipping 2017
        Mp = 0.41*u.Mjup
    Ms = s['st_mass'] * u.Msun
    a = s['pl_orbsmax'] * u.AU
    Rs = s['st_rad'] * u.Rsun
    Rhill = (a * (Mp / 3 / Ms)**(1./3)).to('Rsun')
    Ahill = Rhill*2*Rs*2
    depth = Ahill/(2*pi*Rs**2)
    depth = depth.to('')
    depth[depth > 1] = 1
    return depth


path_etc = reffolder / 'etc_template_1700_1e-13.csv'
etc_ref = table.Table.read(path_etc)
w_etc = etc_ref['wavelength'] * u.AA # EF = earth frame
we_etc = utils.mids2edges(w_etc.value, simple=True) * u.AA
v_etc = w2v(w_etc.value)
pieces = path_etc.name.split('_')
flux_ref = float(pieces[-1][:-4])
expt_ref = float(pieces[-2])
etc_ref['flux2cps'] = etc_ref['target_counts'] / expt_ref / flux_ref
etc_ref['bkgnd_cps'] = (etc_ref['total_counts'] - etc_ref['target_counts'])/expt_ref
def sim_g140m_obs(f, expt):
    fpix = utils.intergolate(we_etc, wgrid_std, f)
    src = fpix * etc_ref['flux2cps'] * expt
    bkgnd = etc_ref['bkgnd_cps'] * expt
    total = bkgnd + src
    with np.errstate(divide='ignore', invalid='ignore'):
        err_counts = np.sqrt(total)
        err_flux = err_counts / src * fpix
        err_flux[~np.isfinite(err_flux)] = 0
    return fpix, err_flux


def doppler_shift(w, velocity):
    return (1 + velocity/const.c)*w


def voigt_xsection(w, w0, f, gamma, T, mass, b=None):
    """
    Compute the absorption cross section using hte voigt profile for a line.

    Parameters
    ----------
    w : astropy quantity array or scalar
        Scalar or vector of wavelengths at which to compute cross section.
    w0: quanitity
        Line center wavelength.
    f: scalar
        Oscillator strength.
    gamma: quantity
        Sum of transition rates (A values) out of upper and lower states. Just Aul for a resonance line where only
        spontaneous decay is an issue.
    T: quantity
        Temperature of the gas. Can be None if you provide a b value instead.
    mass: quantity
        molecular mass of the gas
    b : quantity
        Doppler b value (in velocity units) of the line
    Returns
    -------
    x : quantity
        Cross section of the line at each w.
    """

    nu = const.c / w
    nu0 = const.c / w0
    if T is None:
        sigma_dopp = b/const.c*nu0/np.sqrt(2)
    else:
        sigma_dopp = np.sqrt(const.k_B*T/mass/const.c**2) * nu0
    dnu = nu - nu0
    gauss_sigma = sigma_dopp.to(u.Hz).value
    lorentz_FWHM = (gamma/2/np.pi).to(u.Hz).value
    phi = voigt(dnu.to(u.Hz).value, gauss_sigma, lorentz_FWHM) * u.s
    x = np.pi*const.e.esu**2/const.m_e/const.c * f * phi
    return x.to('cm2')


def voigt(x, gauss_sigma, lorentz_FWHM):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    sigma = gauss_sigma
    gamma = lorentz_FWHM/2.0
    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma /np.sqrt(2*np.pi)
# endregion