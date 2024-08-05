"""This is a mosty unedited copy of code I used for original selection.
I have not updated all paths, only what gets used in scripts32 as of 2024-08-05."""

import os
import glob
from math import pi, sin, cos
import regex as re
from functools import reduce
from copy import copy
pjoin = os.path.join

import numpy as np
from scipy.special import wofz
from matplotlib import pyplot as plt
from astropy import units as u, constants as const, table, coordinates as coord
from tqdm import tqdm
from astroquery.mast import Observations

import galex_motion
import column_map as cmap
import lyasim
import utilities as utils
import paths

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

# region parameters you might care about
Tism = 10000*u.K
# endregion

# # !!!!!!! UPDATE !!!!!!
# catalogs_date = '2024-02'
# folder = '/Users/parke/Google Drive/Research/Proposing/lya transit survey/target search/2024-02 search'
# print('')
# print('HEY! Maker using catalogs from {} folder. Update this parameter in the code if necessary.'.format(catalogs_date))
# print('')


# region utilites

def catpath(name):
    pieces = name.split('.')
    root = ''.join(pieces[:-1])
    suffix = '.' + pieces[-1]
    name_w_date = root + '_' + catalogs_date + suffix
    return pjoin(folder, name_w_date)


def apply_to_quantities(fun, output_unit, *args, **kwargs):
    for i, arg in enumerate(args):
        if isinstance(arg, u.Quantity):
            args[i] = arg.value
        result = fun(*args, **kwargs)
        return result * output_unit


def add_col(cat, name, unit):
    if name in cat.colnames:
        return
    col = table.MaskedColumn(name=name, length=len(cat),
                             dtype='f8')
    col.mask = True
    col.unit = unit
    cat.add_column(col)


def add_src_col(cat, basename):
    srcname = basename + 'src'
    if srcname in cat.colnames:
        return
    src = table.MaskedColumn(name=srcname, length=len(cat),
                             dtype='object')
    src.mask = True
    cat.add_column(src)


def is_pos_real(cat, col):
    if hasattr(cat[col], 'mask'):
        return ~cat[col].mask & (cat[col].filled(-999) > 0)
    else:
        return cat[col] > 0


def Lya_from_galex_schneider19(mag, dist, band='nuv'):
    # mag to flux density from  https://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html
    # Weff from SVO filter service http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=GALEX&asttype=
    dfac = (dist*u.pc/u.AU)**2
    dfac = dfac.to('').value
    if band == 'nuv':
        f = 10**((mag - 20.08)/-2.5)*2.06e-16
        F = f * 768.31
        F1AU = F*dfac
        logFlya = 0.7 * np.log10(F1AU) + 0.193
    elif band == 'fuv':
        f = 10**((mag - 18.82)/-2.5)*1.4e-15
        F = f * 265.57
        F1AU = F*dfac
        logFlya = 0.742 * np.log10(F1AU) + 0.945
    else:
        raise ValueError
    return 10**logFlya


def Lya_from_Teff_linsky13(Teff, Prot):
    fast = Prot < 10
    normal = (Prot >= 10) & (Prot < 25)
    slow = Prot >= 25

    logFlya = np.zeros(len(Teff))

    logFlya[fast] = 0.37688 + 0.0002061*Teff[fast]
    logFlya[normal] = 0.48243 + 0.0001632*Teff[normal]
    logFlya[slow] = -1.5963 + 0.0004732*Teff[slow]

    return 10**logFlya


def EUV_Linsky14(Flya_at_1AU, Teff):
    if type(Teff) is not np.ndarray:
        Teff = np.array([Teff])
        Flya_at_1AU = np.array([Flya_at_1AU])
    logFEUV_bands = np.zeros((9, len(Teff)))
    Ms = Teff < 4000
    logFlya = np.log10(Flya_at_1AU)
    logFEUV_bands[0, Ms] = -0.491 + logFlya[Ms]
    logFEUV_bands[1, Ms] = -0.548 + logFlya[Ms]
    logFEUV_bands[2, Ms] = -0.602 + logFlya[Ms]
    logFEUV_bands[0, ~Ms] = -1.357 + 0.344*logFlya[~Ms] + logFlya[~Ms]
    logFEUV_bands[1, ~Ms] = -1.300 + 0.309*logFlya[~Ms] + logFlya[~Ms]
    logFEUV_bands[2, ~Ms] = -0.882 + logFlya[~Ms]
    logFEUV_bands[3, :] = -2.294+0.258*logFlya + logFlya
    logFEUV_bands[4, :] = -2.098+0.572*logFlya + logFlya
    logFEUV_bands[5, :] = -1.920+0.240*logFlya + logFlya
    logFEUV_bands[6, :] = -1.894+0.518*logFlya + logFlya
    logFEUV_bands[7, :] = -1.811+0.764*logFlya + logFlya
    logFEUV_bands[8, :] = -1.025+ logFlya
    FEUV_bands = 10**logFEUV_bands
    FEUV = np.sum(FEUV_bands, 0)
    return FEUV


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


def sky_distance(ra0, dec0, ra1, dec1):
    """
    Compute arc-distance between two points on sky specified in decimal degrees.
    """
    # convert decimal degrees to radians
    ra0, ra1, dec0, dec1 = map(np.deg2rad, (ra0, ra1, dec0, dec1))

    # haversine formula
    dra = ra1 - ra0
    ddec = dec1 - dec0
    a = np.sin(ddec / 2) ** 2 + np.cos(dec0) * np.cos(dec1) * np.sin(dra / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    return np.rad2deg(c)


def match_by_position(row, cat, planet_or_host='planet'):
    ra, dec = row[['ra', 'dec']]
    ras = cat['ra']
    decs = cat['dec']
    dist = sky_distance(ra, dec, ras, decs)
    close = dist < 0.001
    if planet_or_host == 'host':
        return close
    else:
        P = cat['pl_orbper']
        P = P.filled(0) if hasattr(P, 'filled') else P
        same_period = np.isclose(row['pl_orbper'], P, 0.1)
        return close & same_period


def match_planet_by_name(row, cat):
    name, col = get_id(row)
    if col not in cat.colnames:
        return np.zeros(len(cat), bool)

    result = name == cat[col]
    if hasattr(result, 'filled'):
        return result.filled(False)
    else:
        return result


def get_id(planet):
    cols = ['pl_name', 'toi', 'epic_candname']
    for col in cols:
        if col in planet.colnames:
            if (not getattr(planet[col], 'mask', False)) and planet[col] != '':
                return planet[col], col
    return None


def petigura_gap_lowlim(P):
    logP0 = 1
    logR0 = np.log10(1.84)
    width = 0.1
    m = -0.11
    logR = m*(np.log10(P) - logP0) + logR0 - width/2
    return 10**logR


def ho_gap_lowlim(P, Mstar):
    A, B, C = -0.09, 0.21, 0.35
    width = 0.1
    logP = np.log10(P)
    logM = np.log10(Mstar)
    logR = A*logP + B*logM + C - width/2
    return 10**logR


def rv_transfer(destination_table, source_table):
    """Convenience function for if I want to put RVs from one table into
    another -- like if I collected them by hand for a subset and then want to
    transfer those back into the master search table. """
    cat = destination_table
    for pl in tqdm(source_table):
        printname, _ = get_id(pl)
        matches = match_planet_by_name(pl, cat)
        if sum(matches) == 0:
            matches = match_by_position(pl, cat, 'planet')
        n = sum(matches)
        if n == 0:
            print("No match for {}.".format(printname))
            continue
        if n > 1:
            print("Multiple matches for {}.".format(printname))
            raise ValueError
        i = np.nonzero(matches)[0][0]
        newval = pl['st_radv']
        oldval = cat['st_radv'][i]
        old_is_bad = (not np.isfinite(oldval)) or getattr(oldval, 'mask', False)
        new_is_good = np.isfinite(newval) and (not getattr(newval, 'mask', False))
        if old_is_bad and new_is_good:
            print('Replacing rv for {}, {:.1f} --> {:.1f}'.format(printname, oldval, newval))
            cat['st_radv'][i] = newval


def lya_transit_depth_owen(system_or_systems):
    hnu = 20*u.eV
    photoion_xsection = 6.3e-18*u.cm**2 * (hnu/13.6/u.eV)**-3
    sigma_lya = 2.19250184165e-15*u.cm**2
    eps = 0.1
    vr = 15 * u.km/u.s

    s = system_or_systems
    if isinstance(system_or_systems, table.Row):
        Mp = s['pl_bmasse'] * u.Mearth
        Rp = s['pl_rade'] * u.Rearth
        Ms = s['st_mass'] * u.Msun
        a = s['pl_orbsmax'] * u.AU
        Rs = s['st_rad'] * u.Rsun
        Feuv = s['FEUV_at_planet'] * u.Unit('erg s-1 cm-2')
    else:
        Mp = s['pl_bmasse'].quantity.value * u.Mearth
        Rp = s['pl_rade'].quantity.value * u.Rearth
        Ms = s['st_mass'].quantity.value * u.Msun
        a = s['pl_orbsmax'].quantity.value * u.AU
        Rs = s['st_rad'].quantity.value * u.Rsun
        Feuv = s['FEUV_at_planet'].quantity.value * u.Unit('erg s-1 cm-2')

    Rhill = (a * (Mp / 3 / Ms)**(1./3)).to('Rsun')
    Gamma = Feuv / hnu * photoion_xsection
    Mdot = eps * np.pi * Rp**3 * Feuv / const.G / Mp
    Mdot = Mdot.to('g s-1')
    ncyl = Mdot / (np.pi * Rhill**2 * u.u * vr)
    length = vr / Gamma * np.log(2 * Rhill * sigma_lya * ncyl)
    depth = 4 * Rhill * 2 * length/(pi*Rs**2)
    depth = depth.to('')
    depth[depth > 1] = 1
    return depth


def depth_extended_hill_transit(system_or_systems):
    s = system_or_systems
    if isinstance(system_or_systems, table.Row):
        Mp = s['pl_bmasse'] * u.Mearth
        Ms = s['st_mass'] * u.Msun
        a = s['pl_orbsmax'] * u.AU
        Rs = s['st_rad'] * u.Rsun
    else:
        Mp = s['pl_bmasse'].quantity
        Ms = s['st_mass'].quantity
        a = s['pl_orbsmax'].quantity
        Rs = s['st_rad'].quantity
    Rhill = (a * (Mp / 3 / Ms)**(1./3)).to('Rsun')
    Ahill = Rhill*2*Rs*2
    depth = Ahill/(2*pi*Rs**2)
    depth = depth.to('')
    depth[depth > 1] = 1
    return depth


def unique_star_indices(cat):
    ra, i = np.unique(cat['ra'], return_index=True)
    dec = np.unique(cat['dec'])
    assert np.all(np.sort(cat[i]['dec']) == dec)
    return sorted(i)


def unique_stars(cat):
    return cat[unique_star_indices(cat)]
# endregion


# region catalog processing flow
def make_starting_exocat(include=('composite', 'K2', 'TOI')):
    """
    For most recent files, go retrieve the confirmed planets excluding TESS and K2
    discoveries plus TESS candidates and K2 candidates.

    Returns
    -------

    """
    units_map = {'R_Sun':'solRad',
                 'R_Earth': 'earthRad',
                 'BJD': 'day',
                 'days': 'day',
                 'degrees': 'deg'}

    if 'composite' not in include:
        raise ValueError('Need to at least start with the composite table.')

    table_prefix_map = dict(composite='PSCompPars', K2='k2', TOI='TOI')

    cats = {}
    for key, prefix in table_prefix_map.items():
        if key not in include:
            continue
        files = glob.glob(pjoin(folder, prefix + '*.*.tbl'))
        file = max(files)
        cat = table.Table.read(file, format='ascii.ipac')
        cat = table.Table(cat, masked=True, copy=False)

        # make column names consistent
        cmap.homogenize_columns(cat)

        # make units consistent
        for name in cat.colnames:
            unit = str(cat[name].unit)
            if unit in units_map:
                cat[name].unit = u.Unit(units_map[unit])

        # make trandep and trandur consistent
        if cat['pl_trandep'].unit == 'perc':
            for suff in ('', 'err1', 'err2'):
                cat['pl_trandep' + suff] *= 1e4
                cat['pl_trandep' + suff].unit = 'ppm'
        if cat['pl_trandur'].unit == 'd':
            for suff in ('', 'err1', 'err2'):
                cat['pl_trandur' + suff] *= 24
                cat['pl_trandur' + suff].unit = 'hours'

        if prefix == 'TOI':
            cat['tic_id'] = ['TIC ' + str(id) for id in cat['tic_id']]

        # remove "extra" columns to keep file sizes down (else they take forever to read and write)
        # hope I don't need these later!
        cmap.remove_dispensible_columns(cat)

        # remove known false positives
        if prefix == 'k2':
            bad = cat['disposition'].filled('') == 'FALSE POSITIVE'
        elif prefix == 'TOI':
            bad = cat['tfopwg_disp'].filled('') == 'FP' # false positive
            bad = cat['tfopwg_disp'].filled('') == 'FA' # false alarm
            bad = cat['tfopwg_disp'].filled('') == 'KP' # known planet (should be in other tables)
        else:
            bad = np.zeros(len(cat), bool)
        cat = cat[~bad]

        # remove planets without sky coordinates
        bad = (cat['ra'].mask) | (cat['dec'].mask)
        cat = cat[~bad]

        # remove nontransiting planets
        if prefix == 'PSCompPars':
            nontransiting = (cat['tran_flag'] == False)
            nontransiting = nontransiting.filled(False)
            cat = cat[~nontransiting]

        # remove planets with no orbital period
        cat = cat[~cat['pl_orbper'].mask]

        # set ids
        names = [get_id(row)[0] for row in cat]
        col = table.Column(name='id', length=len(cat), dtype='object')
        cat.add_column(col)
        cat['id'] = names
        assert not np.any(cat['id'] == None)

        cats[key] = cat

    print('Matching and merging catalogs.')
    Nstart = sum([len(tbl) for tbl in cats.values()])

    def match(cat):
        new_cols = list(set(cat.colnames) - set(comp.colnames))
        for name in new_cols:
            col = table.MaskedColumn(name=name, length=len(comp),
                                     dtype=cat[name].dtype, unit=cat[name].unit,
                                     mask=True)
            comp.add_column(col)

        remove = []
        for j, planet in tqdm(list(enumerate(cat))):
            match = match_planet_by_name(planet, comp)
            if not np.any(match):
                match = match_by_position(planet, comp, 'planet')
                if not np.any(match):
                    continue

            if np.sum(match) > 1:
                # these are generally duplicates. Just take the first
                match = np.nonzero(match)[0][0]

            remove.append(j)
            for name in new_cols:
                comp[name][match] = cat[name][j]

        cat.remove_rows(remove)

    comp = cats['composite']
    if 'K2' in include:
        k2 = cats['K2']
        print('Matching K2 candidates.')
        match(k2)
        comp = table.vstack((comp, k2))
        comp = table.Table(comp, masked=True) # bug in asstropy that vstack unmasks?
        # don't stack all three at once or else some K2 and TESS duplicates might not be matched

    if 'TOI' in include:
        tess = cats['TOI']
        print('Matching TESS candidates.')
        match(tess)
        comp = table.vstack((comp, tess))

    print('Started with {} planets (after filtering out bad apples).'.format(Nstart))
    print('Combining catalogs resulted in {}.'.format(len(comp)))

    # remove duplicates from the composite table
    print('Removing duplicates from the merged catalog.')
    remove = []
    duplicates = []
    icheck = list(range(len(comp)))
    while len(icheck):
        i = icheck[0]
        print('\r{}/{}'.format(i, len(comp)), end="", flush=True)
        matches = match_by_position(comp[i], comp, 'planet')
        j, = np.nonzero(matches)
        if len(j) > 1:
            duplicates.extend(j)
        remove.extend(j[1:])
        [icheck.remove(jj) for jj in j]
    print('The combined catalog had these duplicates:')
    comp[duplicates]['id ra dec pl_orbper pl_rade'.split()].pprint(-1)
    comp.remove_rows(remove)
    print('Removing duplicates resulted in {}'.format(len(comp)))

    print('Saving the new catalog. This will take a few minutes... be patient.')
    comp.write(catpath('0_merged_starting_catalog.ecsv'))
    radec = comp[['ra', 'dec']]
    rapath = catpath('merged_starting_catalog_radec.csv')
    # need the date in there cuz x-match service doesn't like duplicate filenames
    np.savetxt(rapath, radec, delimiter=',', header='ra,dec')

    return comp


def add_simbad(cat, simbad):
    """
    Match the ra and dec saved as radec_date.csv when you ran make_starting_exocat to SIMBAD via
    CDS X-match. Put those two tables into this. Consider using a match distance larger than the default (say 60").

    """

    combo = cat.copy()
    for name in simbad.colnames:
        if name not in combo.colnames:
            col = table.MaskedColumn(name=name, length=len(combo),
                                     dtype=simbad[name].dtype)
            col.mask = True
            combo.add_column(col)
    col = table.MaskedColumn(name='multimatch', length=len(combo), dtype=bool)
    combo.add_column(col)

    for k, row in tqdm(list(enumerate(combo))):
        i, = np.nonzero(np.isclose(row['ra'], simbad['ra']))
        if len(i) > 1:
            d = simbad['d_arcsec'][i]
            j = np.argmin(d)
            match = i[j]
            multi = True
        elif len(i) == 1:
            match = i[0]
            multi = False
        else:
            continue

        for name in simbad.colnames:
            combo[name][k] = simbad[name][match]
        combo['multimatch'][k] = multi

    # move rvs over
    add_src_col(combo, 'st_radv')
    move = combo['st_radv'].mask & ~combo['radvel'].mask
    already_good = ~combo['st_radv'].mask
    combo['st_radvsrc'][already_good] = 'NASA Exoplanet Archive'
    combo['st_radv'][move] = combo['radvel'][move]
    combo['st_radvsrc'][move] = 'SIMBAD'

    # move pms over
    add_src_col(combo, 'sy_pmra')
    move = combo['sy_pmra'].mask & ~combo['pmra'].mask
    already_good = ~combo['sy_pmra'].mask
    combo['sy_pmrasrc'][already_good] = 'NASA Exoplanet Archive'
    combo['sy_pmra'][move] = combo['pmra'][move]
    combo['sy_pmrasrc'][move] = 'SIMBAD'

    # move plx over
    add_src_col(combo, 'sy_plx')
    move = combo['sy_plx'].mask & ~combo['plx'].mask
    already_good = ~combo['sy_plx'].mask
    combo['sy_plxsrc'][already_good] = 'NASA Exoplanet Archive'
    combo['sy_plx'][move] = combo['plx'][move]
    combo['sy_plxsrc'][move] = 'SIMBAD'

    # compute and move plx based distance
    dist = 1000/combo['plx']
    add_src_col(combo, 'sy_dist')
    move = combo['sy_dist'].mask & ~combo['plx'].mask
    already_good = ~combo['sy_dist'].mask
    combo['sy_distsrc'][already_good] = 'NASA Exoplanet Archive'
    combo['sy_dist'][move] = dist[move]
    combo['sy_distsrc'][move] = 'SIMBAD'

    # move Jmag over (for TSM values calculated later, although I think exocat probs has them all)
    move = combo['sy_jmag'].mask & ~combo['J'].mask
    combo['sy_jmag'][move] = combo['J'][move]

    add_src_col(combo, 'sy_pmdec')
    move = combo['sy_pmdec'].mask & ~combo['pmdec'].mask
    already_good = ~combo['sy_pmdec'].mask
    combo['sy_pmdecsrc'][already_good] = 'NASA Exoplanet Archive'
    combo['sy_pmdec'][move] = combo['pmdec'][move]
    combo['sy_pmdecsrc'][move] = 'SIMBAD'

    cmap.remove_dispensible_columns(combo)

    path = catpath('1_planets_w_simbad.ecsv')
    combo.write(path)
    return combo


def add_sam_rvs(cat, samrvs):
    """
    You can probs read in samrvs just by setting format='ascii' in
    table.read. If you have rvs from both chiron and tres, just stack those
    tables before feeding to this.

    Be sure to check if any appealing rv targets have large std devs on rv in
    sam's table, as this could indicate a binary.

    Parameters
    ----------
    cat
    samrvs

    Returns
    -------

    """
    for row in tqdm(samrvs):
        toi = row['TOI']
        match = cat['toipfx'].filled('0').astype('int') == toi
        cat['st_radv'][match] = row['<V>']
        cat['st_radvsrc'][match] = 'TRES from sam'

    print("Be sure to check if any appealing rv targets have large std devs on rv in sam's table, as this could indicate a binary.")

    path = catpath('2_planets_w_sam_rvs.ecsv')
    cat.write(path)

    return cat


def merge_past_rvs(newcat, oldcat):
    # add rvs to the norvcat from the previous function and this will merge them into the combo catalog
    # I used transnr > 1.5 as a cutoff in the past
    # then you can rerun the previous function to be sure you have all good systems

    new = newcat.copy()
    rv_transfer(new, oldcat) # this happens in place
    path = catpath('3_planets_w_old_by_hand_rvs.ecsv')
    new.write(path)

    return new


def groom_evapgap_galex(evapgap_galex):
    cat = evapgap_galex
    new = cat.copy()

    nuv_cps = cat['galex_nuv']
    nuv_cpserr = cat['galex_nuverr1']
    nuv_mag = galex_motion.galex_cps2mag(nuv_cps, 'nuv')
    nuv_magerr = galex_motion.cps2magerr(nuv_cps, nuv_cpserr)
    new['galex_nuv'] = nuv_mag
    new['galex_nuverr1'] = nuv_magerr
    new['galex_nuverr2'] = nuv_magerr
    new['galex_nuvlim'] = -cat['galex_nuvlim']

    fuv_cps = cat['galex_fuv']
    fuv_cpserr = cat['galex_fuverr1']
    fuv_mag = galex_motion.galex_cps2mag(fuv_cps, 'fuv')
    fuv_magerr = galex_motion.cps2magerr(fuv_cps, fuv_cpserr)
    new['galex_fuv'] = fuv_mag
    new['galex_fuverr1'] = fuv_magerr
    new['galex_fuverr2'] = fuv_magerr
    new['galex_fuvlim'] = -cat['galex_fuvlim']

    return new


def add_galex_cols_inplace(cat):
    add_col(cat, 'sy_fuvmag', None)
    add_col(cat, 'sy_fuvmagerr1', None)
    add_col(cat, 'sy_fuvmagerr2', None)
    add_col(cat, 'sy_fuvmaglim', None)

    add_col(cat, 'sy_nuvmag', None)
    add_col(cat, 'sy_nuvmagerr1', None)
    add_col(cat, 'sy_nuvmagerr2', None)
    add_col(cat, 'sy_nuvmaglim', None)


def merge_past_galex(destination, source, colname_root='galex_', exception_handling='return'):
    cat = destination.copy()

    add_galex_cols_inplace(cat)

    def transfer(i, j, destcol, srccol):
        newval = source[srccol][j]
        oldval = cat[destcol][i]
        old_is_bad = (not np.isfinite(oldval)) or getattr(oldval, 'mask', False)
        new_is_good = np.isfinite(newval) and (not getattr(newval, 'mask', False))
        if old_is_bad and new_is_good:
            cat[destcol][i] = newval
            cat[destcol+'err1'][i] = source[srccol+'err1'][j]
            cat[destcol+'err2'][i] = source[srccol+'err2'][j]
            cat[destcol+'lim'][i] = source[srccol+'lim'][j]

    try:
        for i, pl in tqdm(list(enumerate(cat))):
            matches = match_planet_by_name(pl, source)
            if sum(matches) == 0:
                matches = match_by_position(pl, source, 'planet')
            n = sum(matches)
            if n == 0:
                continue
            j = np.nonzero(matches)[0][0] # note this could get buggy if there are multiple matches
            transfer(i, j, 'sy_fuvmag', colname_root+'fuv')
            transfer(i, j, 'sy_nuvmag', colname_root+'nuv')
    except:
        if exception_handling == 'return':
            return cat
        else:
            raise

    return cat


def add_galex(cat, exception_handling="return", timeout=600., write=False):
    # probably because this is slow on the server side you will want to divide up the table and have multiple threads each doing a piece
    # in which case you will want write to be false

    # newcat = cat.copy()  # with astropy 4.2.1 this is causing errors later when the table is saved
    raise(ValueError('You need to update this to use proper motions from simbad because archive ones use weird epoch.'))
    newcat = cat

    add_galex_cols_inplace(newcat)

    try:
        for i, row in tqdm(list(enumerate(newcat))):
            if ~newcat['sy_nuvmag'].mask[i]:
                continue

            ra = row['ra']
            dec = row['dec']
            if newcat['sy_pmra'].mask[i]:
                pm_ra = 0
                pm_dec = 0
            else:
                pm_ra = row['sy_pmra']
                pm_dec = row['sy_pmdec']

            result = galex_motion.extract_and_coadd(ra, dec, pm_ra, pm_dec,
                                                    match_radius=16./3600,
                                                    query_timeout=timeout)

            (nuv, nerr), (fuv, ferr) = result
            if nuv < 0:
                newcat['sy_nuvmag'][i] = galex_cps2mag(nerr, 'nuv')
                newcat['sy_nuvmaglim'][i] = -1
            else:
                newcat['sy_nuvmag'][i] = galex_cps2mag(nuv, 'nuv')
                magerr = galex_motion.cps2magerr(nuv, nerr)
                newcat['sy_nuvmagerr1'][i] = magerr
                newcat['sy_nuvmagerr2'][i] = magerr

            if fuv < 0:
                newcat['sy_fuvmag'][i] = galex_cps2mag(ferr, 'nuv')
                newcat['sy_fuvmaglim'][i] = -1
            else:
                newcat['sy_fuvmag'][i] = galex_cps2mag(fuv, 'fuv')
                magerr = galex_motion.cps2magerr(fuv, ferr)
                newcat['sy_fuvmagerr1'][i] = magerr
                newcat['sy_fuvmagerr2'][i] = magerr
    except:
        if exception_handling == 'return':
            return newcat
        else:
            raise

    def set_nan_mask(col):
        nans = np.isnan(newcat[col])
        newcat[col].mask[nans] = True
    map(set_nan_mask, ('sy_nuvmag', 'sy_nuvmagerr1', 'sy_nuvmagerr2',
                       'sy_fuvmag', 'sy_fuvmagerr1', 'sy_fuvmagerr2'))

    if write:
        path = catpath('1.5_planets_w_galex.ecsv')
        newcat.write(path)

    return newcat


def fill_in_parameters(cat, save=True):

    # tools for replacing empties with estimates
    isgood = lambda col: is_pos_real(cat, col)
    replace = lambda a,b: isgood(a) & ~isgood(b)

    # now fill in stellar parameters (mass, Teff, radius)
    st_cols = ('st_rad', 'st_teff', 'st_mass', 'sy_dist')
    for col in st_cols:
        add_src_col(cat, col)
        cat[col + 'src'][isgood(col)] = 'exocat'
    ms = table.Table.read(paths.reffolder / 'main_sequence.dat', format='ascii.csv')
    ms.sort('Teff')

    ## distance from parallax
    rep = replace('sy_plx', 'sy_dist')
    d = 1000 / cat['sy_plx'] * u.pc
    cat['sy_dist'][rep] = d[rep]
    cat['sy_distsrc'][rep] = 'simbad'

    # using Teff
    R = np.interp(cat['st_teff'], ms['Teff'], ms['R'])
    rep = replace('st_teff', 'st_rad')
    cat['st_rad'][rep] = R[rep]
    cat['st_radsrc'][rep] = 'interp from st_teff'

    M = np.interp(cat['st_teff'], ms['Teff'], ms['M'])
    rep = replace('st_teff', 'st_mass')
    cat['st_mass'][rep] = M[rep]
    cat['st_masssrc'][rep] = 'interp from st_teff'

    # using R
    # FIXME update to this table https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
    T = np.interp(cat['st_rad'], ms['R'], ms['Teff'])
    rep = replace('st_rad', 'st_teff')
    cat['st_teff'][rep] = T[rep]
    cat['st_teffsrc'][rep] = 'interp from st_rad'

    M = np.interp(cat['st_rad'], ms['R'], ms['M'])
    rep = replace('st_rad', 'st_mass')
    cat['st_mass'][rep] = M[rep]
    cat['st_masssrc'][rep] = 'interp from st_rad'

    # now fill in planet mass, eqt, and a
    pl_cols = ('pl_orbsmax', 'pl_eqt', 'pl_bmasse')
    for col in pl_cols:
        add_src_col(cat, col)
        cat[col + 'src'][isgood(col)] = 'exocat'

    ## a using stellar Teff, rad, and planet Teq
    rep = isgood('st_teff') & isgood('st_rad') & isgood('pl_eqt') & ~isgood('pl_orbsmax')
    Tstar = cat['st_teff'].quantity.value
    Tplanet = cat['pl_eqt'].quantity.value
    Rstar = cat['st_rad'].quantity.value * u.Rsun
    assert np.all(Tplanet[rep] > 0)
    with np.errstate(divide='ignore'):
        a = (Tstar/Tplanet)**2 * Rstar
    a = a.to('AU')
    cat['pl_orbsmax'][rep] = a[rep]
    cat['pl_orbsmaxsrc'][rep] = 'from pl_eqt, st_teff, and st_rad'

    ## a using period and stellar mass
    rep = isgood('st_mass') & isgood('pl_orbper') & ~isgood('pl_orbsmax')
    P = cat['pl_orbper'].quantity.value * u.d
    M = cat['st_mass'].quantity.value * u.Msun
    G = const.G
    a = (G*M/4/pi**2 * P**2)**(1/3)
    a = a.to('AU')
    cat['pl_orbsmax'][rep] = a[rep]
    cat['pl_orbsmaxsrc'][rep] = 'from st_mass and pl_orbper'

    ## Teq using Lstar and a
    rep = isgood('st_lum') & isgood('pl_orbsmax') & ~isgood('pl_eqt')
    L = 10 ** cat['st_lum'].quantity.value * u.Lsun
    a = cat['pl_orbsmax'].quantity.value * u.AU
    sb = const.sigma_sb
    Teq = (L/16/pi/sb/a**2)**(1/4)
    Teq = Teq.to('K')
    cat['pl_eqt'][rep] = Teq[rep]
    cat['pl_eqtsrc'][rep] = 'from st_lum and pl_orbsmax'

    ## Teq using Teff and Rstar and a
    rep = (isgood('st_teff') & isgood('st_rad') & isgood('pl_orbsmax')
           & ~isgood('pl_eqt'))
    Teff = cat['st_teff'].quantity.value * u.K
    Rstar = cat['st_rad'].quantity.value * u.Rsun
    a = cat['pl_orbsmax'].quantity.value * u.AU
    Teq = np.sqrt(Rstar/2/a) * Teff
    Teq = Teq.to('K')
    cat['pl_eqt'][rep] = Teq[rep]
    cat['pl_eqtsrc'][rep] = 'from st_teff, st_rad, and pl_orbsmax'

    ## Mplanet from Rplanet
    # note that whenever there is a massj there is a masse -- I checked
    rep = isgood('pl_rade') & ~isgood('pl_bmasse')
    Rp = cat['pl_rade'].quantity.value
    small = Rp < 1.23
    medium = ~small & (Rp < 14.26)
    large = Rp >= 14.26
    M = np.zeros_like(Rp)
    M[small] = 0.9718 * Rp[small]**3.58
    M[medium] = 1.436 * Rp[medium]**1.7
    rep[large] = False # after the breakpoint, increases in mass hardly change radius
    cat['pl_bmasse'][rep] = M[rep]
    cat['pl_bmassesrc'][rep] = 'from pl_rade per Chen & Kipping 2017'

    pairs = (('pl_eqt', '%.0f'),
             ('pl_orbsmax', '%.3f'),
             ('pl_rade', '%.3f'),
             ('pl_bmasse', '%.3f'),
             ('pl_orbper', '%.3f'),
             ('st_rad', '%.3f'),
             ('st_mass', '%.3f'),
             ('sy_dist', '%.3f'),
             ('st_teff', '%.0f'))
    for name, fmt in pairs:
        cat[name].format = fmt

    if save:
        path = catpath('4_planets_w_parmeters_filled.ecsv')
        cat.write(path)

    return cat


def fill_in_Lya_EUV_TSM_Q(cat, save=True):
    isgood = lambda col: is_pos_real(cat, col)

    add_col(cat, 'Flya_at_1AU', u.Unit('erg s-1 cm-2'))
    filled = np.zeros(len(cat), bool)

    # Lya based on GALEX NUV
    if 'NUV' in cat.colnames:
        K_or_M = cat['st_teff'] < 4500. # because schneider GALEX-Lya relationship is for late Ks and Ms only
        rep = isgood('NUV') & isgood('sy_dist') & K_or_M & ~filled
        filled = filled | rep
        print('Using GALEX NUV per Schneider+ 2019 to predict Lya for {} rows.'.format(sum(rep)))
        Flya_au = Lya_from_galex_schneider19(cat['NUV'], cat['sy_dist'], 'nuv')
        Flya_au *= u.Unit('erg s-1 cm-2')
        cat['Flya_at_1AU'][rep] = Flya_au[rep]

    # Lya based on Teff from Linsky+ 2013
    rep = isgood('st_teff') & isgood('st_rotp') & ~filled
    filled = filled | rep
    print('Using Teff per Linsky+ 2013 to predict Lya for {} rows.'.format(sum(rep)))
    Flya_au = Lya_from_Teff_linsky13(cat['st_teff'], cat['st_rotp'])
    Flya_au *= u.Unit('erg s-1 cm-2')
    cat['Flya_at_1AU'][rep] = Flya_au[rep]

    # Lya based on Teff from Schneider
    rep = isgood('st_teff') & ~filled
    print('Using Teff per Schneider+ 2019 to predict Lya for {} rows.'.format(sum(rep)))
    Teff = cat['st_teff'].quantity.value
    logF = 6.754e-4 * Teff - 2.639
    Flya_au = 10**logF * u.Unit('erg s-1 cm-2')
    cat['Flya_at_1AU'][rep] = Flya_au[rep]

    Flya_au = cat['Flya_at_1AU'].quantity

    add_col(cat, 'Flya_at_earth', u.Unit('erg s-1 cm-2'))
    rep = isgood('Flya_at_1AU') & isgood('sy_dist')
    d = cat['sy_dist'].quantity.value * u.pc
    assert np.all(d[rep] > 0)
    with np.errstate(divide='ignore'):
        Flya_earth = Flya_au * (u.AU/d)**2
    Flya_earth = Flya_earth.to('erg s-1 cm-2')
    cat['Flya_at_earth'][rep] = Flya_earth[rep]

    add_col(cat, 'Flya_at_planet', u.Unit('erg s-1 cm-2'))
    rep = isgood('Flya_at_1AU') & isgood('pl_orbsmax')
    a = cat['pl_orbsmax'].quantity.value * u.AU
    Flya_planet = Flya_au * (u.AU / a) ** 2
    Flya_planet = Flya_planet.to('erg s-1 cm-2')
    cat['Flya_at_planet'][rep] = Flya_planet[rep]

    print('Estimating EUV from Lya per Linsky+ 2014')

    add_col(cat, 'FEUV_at_1AU', u.Unit('erg s-1 cm-2'))
    FEUV = EUV_Linsky14(Flya_au.value, Teff) * u.Unit('erg s-1 cm-2')
    cat['FEUV_at_1AU'][rep] = FEUV[rep]

    add_col(cat, 'FEUV_at_earth', u.Unit('erg s-1 cm-2'))
    rep = isgood('FEUV_at_1AU') & isgood('sy_dist')
    d = cat['sy_dist'].quantity.value * u.pc
    assert np.all(d[rep] > 0)
    with np.errstate(divide='ignore'):
        FEUV_earth = FEUV * (u.AU/d)**2
    FEUV_earth = FEUV_earth.to('erg s-1 cm-2')
    cat['FEUV_at_earth'][rep] = FEUV_earth[rep]

    add_col(cat, 'FEUV_at_planet', u.Unit('erg s-1 cm-2'))
    rep = isgood('FEUV_at_1AU') & isgood('pl_orbsmax')
    a = cat['pl_orbsmax'].quantity.value * u.AU
    FEUV_planet = FEUV * (u.AU / a) ** 2
    FEUV_planet = FEUV_planet.to('erg s-1 cm-2')
    cat['FEUV_at_planet'][rep] = FEUV_planet[rep]

    # region TSM
    add_col(cat, 'TSM', u.Unit(''))
    rep = (isgood('pl_rade') & isgood('pl_eqt') & isgood('pl_bmasse')
           & isgood('st_rad') & isgood('sy_jmag'))
    Rp = cat['pl_rade'].quantity.value
    Teq = cat['pl_eqt'].quantity.value
    Mp = cat['pl_bmasse'].quantity.value
    Rs = cat['st_rad'].quantity.value
    J = cat['sy_jmag'].quantity.value
    R1 = (Rp < 1.5)
    R2 = (Rp >= 1.5) & (Rp < 2.)
    R3 = (Rp >= 2) & (Rp < 2.75)
    R4 = (Rp >= 2.75) & (Rp < 4)
    R5 = (Rp >= 4) & (Rp < 10.)
    TSM_raw = Rp**3*Teq/Mp/Rs**2 * 10**(-J/5)
    cat['TSM'][R1 & rep] = TSM_raw[R1 & rep] * 0.19
    cat['TSM'][R2 & rep] = TSM_raw[R2 & rep] * 1.26 * 2.3/18
    cat['TSM'][R3 & rep] = TSM_raw[R3 & rep] * 1.26
    cat['TSM'][R4 & rep] = TSM_raw[R4 & rep] * 1.28
    cat['TSM'][R5 & rep] = TSM_raw[R5 & rep] * 1.15
    # endregion

    # region ohmic heating
    norm = 10**-1.4 * u.W/u.m**2
    # norm is set to match the middle of the sig = 0.1 and 10 cases for di = 1000 u.km in Fig 4 of Cohen+ 2024
    # for the reference values given below solar case, 0.1 AU
    a_ref = 0.1
    dBdt_ref = 0.1
    P_ref = 0.1**(3/2)
    dBdt = dBdt_ref * (cat['pl_orbsmax']/a_ref)**-3 * (cat['pl_orbper']/P_ref)**-1
    # the above assumes all stars have similar mag field and it scales as dipole
    # (supported somewhat by sun and trappist cases in Cohen+ 2024)
    Q = norm * dBdt.value**2
    cat['joule_heating'] = Q.to('erg s-1 cm-2')
    # endregion

    pairs = (('Flya_at_1AU', '%.2e'),
             ('Flya_at_earth', '%.2e'),
             ('Flya_at_planet', '%.2f'),
             ('FEUV_at_1AU', '%.2e'),
             ('FEUV_at_earth', '%.2e'),
             ('FEUV_at_planet', '%.2f'),
             ('TSM', '%.1f'))
    for name, fmt in pairs:
        cat[name].format = fmt

    if save:
        path = catpath('5_planets_w_Lya_estimates.ecsv')
        cat.write(path)

    return cat


def add_transit_snr(cat):
    # estimate SNR of transit in absorbed ISM profile
    n_H_nom = 0.031 / u.cm ** 3  # mid value
    n_H_opt = 0.01/u.cm**3
    Tism = 10000 * u.K
    wgrid = np.arange(1210., 1220., 0.005) * u.AA

    etc_flux = 1e-14
    etc_src = 105 * 2 / 10000.
    etc_bg = 10.5 * 2 / 10000.

    v_integrate = [-150, 100]

    v_trans, d_trans = np.loadtxt(paths.reffolder / 'generic profile.csv', delimiter=',').T
    d_trans_norm = d_trans / np.max(d_trans)

    cat['rv_star_ism'] = 0.0
    cat['flag_lya_core_exposed'] = False
    pairs = [['nominal', n_H_nom],
             ['optimistic', n_H_opt]]
    for name, n_H in pairs:
        snr_col = 'pl_lya_tnst_snr_{}'.format(name)
        cat[snr_col] = 0.0
        for i in range(0, len(cat)):
            star = cat[i]
            # scale depth to match extended hill sphere transit
            max_depth = depth_extended_hill_transit(star)
            d_trans = max_depth * d_trans_norm
            Nh = n_H*star['sy_dist'] * u.pc
            rv_ism = ism_velocity(star['ra'], star['dec'])
            rv_star = star['st_radv'] * u.km/u.s
            if np.ma.is_masked(rv_star):
                rv_star = rv_ism # assume no offset between star and ISM
            cat['rv_star_ism'][i] = rv_star.value - rv_ism.value
            lya = star['Flya_at_earth'] * u.Unit('erg s-1 cm-2')
            v = (wgrid/wlab_H - 1)*3e5 - rv_star.value

            intrinsic = reversed_lya_profile(wgrid, rv_star, lya)
            transmitted = transmission(wgrid, rv_ism, Nh, Tism)
            observed = intrinsic * transmitted
            flux_core = np.interp(0, v, observed)
            if flux_core.value > 1e-14:
                cat['flag_lya_core_exposed'][i] = True
            transit_depth = np.interp(v, v_trans, d_trans)
            transit = observed * (1 - transit_depth)

            # FIXME need to fix this. it doesn't bin as the etc does!! SNR is way overestimated
            src_cps = observed.value/etc_flux * etc_src
            src_cnts = src_cps* 3 * 2700.
            bg_cnts = etc_bg * 3 * 2700.
            err = np.sqrt(src_cnts + bg_cnts)/src_cnts*observed
            err[np.isnan(err)] = 0

            mask = (v > v_integrate[0]) & (v < v_integrate[1])
            Fout = np.trapz(observed[mask], v[mask])
            Fin = np.trapz(transit[mask], v[mask])
            Eout = np.sqrt(np.trapz(err[mask]**2, v[mask]))
            dF = Fout - Fin
            dE = np.sqrt(2)*Eout
            SNR = dF/dE
            cat[snr_col][i] = SNR

    return cat


def initial_cut(cat):
    mask = cat['pl_lya_tnst_snr_optimistic'] > 3
    return cat[mask]


def flag_jwst_hst_spectra(cat):
    cat['flag_jwst_spectra'] = False
    cat['flag_hst_spectra'] = False
    ra_searched = []
    for i, planet in tqdm(list(enumerate(cat))):
        ra = planet['ra']
        dec = planet['dec']
        if ra not in ra_searched: # no need to duplicate
            n = Observations.query_criteria_count(obs_collection='JWST',
                                                  coordinates=f"{ra} {dec:+f}",
                                                  radius='60s',
                                                  dataproduct_type='spectrum')
            if n > 0:
                cat['flag_jwst_spectra'][i] = True
            n = Observations.query_criteria_count(obs_collection='HST',
                                                  coordinates=f"{ra} {dec:+f}",
                                                  radius='60s',
                                                  dataproduct_type='spectrum',
                                                  wavelength_region=['Optical', 'Infrared'])
            # fixme should probably check abstracts for mention of transit
            if n > 0:
                cat['flag_hst_spectra'][i] = True
        ra_searched.append(ra)


def flag_rv_systems(cat):
    import etta

    # first just look for planets with masses already measured
    cat['flag_goodmass'] = cat['pl_bmassprov'] == 'Mass'
    cat['flag_tfop_rvs'] = False

    # next look for precision RV data in ExoFOP if planet is a TOI
    tics = []
    for i, planet in tqdm(list(enumerate(cat))):
        tic = planet['tic_id'][4:]
        if not np.ma.is_masked(tic):
            if tic in tics:
                continue
            tbl = etta.download_spect(target=tic)
            n = np.sum(tbl['Appropriate to PRV'] == 'Yes')
            planet['flag_tfop_rvs'] = n > 0
            tics.append(tic)


def flag_interesting_systems(cat):
    cat['flag_snr_okay'] = snr_okay = cat['pl_lya_tnst_snr_nominal'] > 3
    cat['flag_snr_xcptn'] = snr_xcptn = cat['pl_lya_tnst_snr_optimistic'] > 3

    # only planets for which hydro escaping H/He atmospheres are realistic
    Rgap = ho_gap_lowlim(cat['pl_orbper'], cat['st_mass'])
    cat['flag_abovegap'] = abovegap = cat['pl_rade'] > Rgap.filled(1.8)
    dR_gap_dex = np.log10(cat['pl_rade']/Rgap.filled(100))
    cat['flag_gap_upper_cusp'] = (dR_gap_dex > 0) & (dR_gap_dex < 0.1)
    cat['flag_dense'] = dense = cat['pl_dens'].filled(0) > 4
    cat['flag_warm'] = warm = cat['pl_eqt'].filled(500) > 300.
    cat['flag_young'] = young = cat['st_age'].filled(5) < 1
    # lower xuv threshold to catch stars where code likely underpredicts xuv flux
    cat['flag_hydro'] = hydro = (cat['FEUV_at_planet'] > 1) | (young & (cat['FEUV_at_planet'] > 0.1))
    # actually, there are only 11 planets with ~hydro,
    # and I think we should keep them to test when escape turns off,
    # so no longer selecting for hydro
    cat['flag_HHe'] = HHe = abovegap.filled(False) & ~dense

    # I'll make these chosen lists and then merge them later
    chosens = [HHe & snr_okay]
    chosens.append(young & HHe & snr_xcptn)

    # multis
    multi = np.zeros(len(cat), bool)
    multi_hydro = np.zeros(len(cat), bool)
    multi_gap = np.zeros(len(cat), bool)
    for planet in cat:
        maskra = cat['ra'] == planet['ra']
        maskdec = cat['dec'] == planet['dec']
        assert np.all(maskra == maskdec)
        m = maskra
        if np.sum(m) > 1:
            multi[m] = True
            hydro_on = HHe[m] & hydro[m]
            hydro_off = HHe[m] & ~hydro[m]
            multi_hydro[m] = (np.sum(hydro_on) > 0) & (np.sum(hydro_off) > 0)
            above_gap = abovegap[m] & hydro[m]
            below_gap = ~abovegap[m] & hydro[m]
            multi_gap[m] = (np.sum(above_gap) > 0) & (np.sum(below_gap) > 0)
    cat['flag_multi'] = multi
    cat['flag_multi_hydro_onoff'] = multi_hydro
    cat['flag_multi_above_below_gap'] = multi_gap
    chosens.append(multi_hydro & snr_xcptn)
    chosens.append(multi_gap & snr_xcptn)

    # water world compatible
    mm, rr = np.loadtxt('../auxilliary data/lp22 m-r water worlds.csv', delimiter=',').T
    pp = mm*u.Mearth/(4/3*np.pi*(rr*u.Rearth)**3)
    pp = pp.to_value('g cm-3')
    # mme, rre = np.loadtxt('../auxilliary data/lp22 m-r earth.csv', delimiter=',').T
    # ppe = mme*u.Mearth/(4/3*np.pi*(rre*u.Rearth)**3)
    # ppe = ppe.to_value('g cm-3')
    # plt.loglog(mme, ppe)
    # plt.loglog(mm, pp)
    # plt.show()
    # goodmass = combo['pl_bmasseerr1']/combo['pl_bmasse'] < 0.25
    # goodrad = combo['pl_radeerr1']/combo['pl_rade'] < 0.08
    p_water = np.interp(cat['pl_bmasse'], mm, pp, left=np.nan, right=np.nan)
    rho_diff_sig = np.abs((cat['pl_dens'] - p_water) / cat['pl_denserr1'])
    water_world = (rho_diff_sig < 1) & (cat['pl_bmasse'] < 30) & (cat['pl_rade'] < 3)
    water_world = water_world.filled(False)
    cat['flag_water_world'] = water_world
    chosens.append(water_world & snr_xcptn)

    # super puffs
    p_sp = 0.3 # based on Libby-Roberts 2020
    super_puff = (cat['pl_dens'].filled(10) < p_sp) & (cat['pl_bmasse'] < 30)
    cat['flag_super_puff'] = super_puff
    chosens.append(super_puff.filled(False) & snr_xcptn)
    # there are none of these in the sample. too bad.


    # HZ planets
    # this might mainly select for young planets
    ai, mi = np.loadtxt('../auxilliary data/kop hz inner.csv', delimiter=',', unpack=True)
    ao, mo = np.loadtxt('../auxilliary data/kop hz outer.csv', delimiter=',', unpack=True)
    ai_p = 10**np.interp(np.log10(cat['st_mass']), np.log10(mi), np.log10(ai))
    ao_p = 10**np.interp(np.log10(cat['st_mass']), np.log10(mo), np.log10(ao))
    cat['flag_hz'] = hz = (cat['pl_orbsmax'] > ai_p) & (cat['pl_orbsmax'] < ao_p)
    chosens.append(hz.filled(False) & snr_xcptn)


    # high TSM
    cat['flag_high_TSM'] = high_TSM = (((cat['pl_rade'] <= 1.5) & (cat['TSM'] > 10))
                                         | ((cat['pl_rade'] > 1.5) & (cat['TSM'] > 90)))
    chosens.append(high_TSM.filled(False) & snr_xcptn)

    # high eccentricity
    # cat['flag_eccentric'] = eccentric = cat['pl_orbeccen'] > 0.1
    # chosens.append(eccentric & snr_xcptn)

    # joule heating
    # I estimated this and it was never above 1% even in extreme cases for planets

    # systems with other observtions/follow up add to chosens
    chosens.append(cat['flag_jwst_spectra'] & snr_xcptn)
    chosens.append(cat['flag_hst_spectra'] & snr_xcptn)
    chosens.append(cat['flag_goodmass'].filled(False) & snr_xcptn)
    chosens.append(cat['flag_tfop_rvs'] & snr_xcptn)

    # merge chosens
    chosen = reduce(np.logical_or, chosens[1:], chosens[0])
    cat['flag_chosen'] = chosen
    return chosens


def flag_duplicates(cat):

    # region read in abstracts
    file = folder + '/abstracts.cat'
    with open(file) as f:
        txt = f.read()
    abstracts = txt.split('------------------------------------------------------------------------------')
    abstracts = abstracts[:-1]

    abdic, fails = {}, []
    for a in abstracts:
        result = re.search('ID: +(\d+)\n', a)
        if result is None:
            fails.append(a)
            continue
        id, = result.groups(0)
        id = int(id)
        abdic[id] = a
    # endregion

    # region read in obesrvations
    file = folder + '/paec_7-present.cat' # be careful not to download the _ss_ catalog as I think it is solar system objects
    obs = table.Table.read(file, format='ascii.fixed_width_two_line', position_line=18, header_start=17)
    # some basic cuts to reduce the table size
    keep = ((np.char.count(obs['config'], 'STIS') > 0)
            | (np.char.count(obs['config'], 'COS') > 0)
            & (np.char.count(obs['mode'], 'ACQ') == 0))
    obs = obs[keep]
    obs['ra'] = coord.Angle(obs['ra'], unit=u.hourangle)
    obs['dec'] = coord.Angle(obs['dec'], unit=u.deg)
    obs_coords = coord.SkyCoord(obs['ra'], obs['dec'])
    # endregion

    # region search for observations
    n = len(cat)
    cat['hst abstracts'] = table.MaskedColumn('', dtype=object, length=n, fill_value='', mask=True)
    cat['n_e140m_obs'] = table.MaskedColumn(0, dtype=int, length=n)
    cat['n_g140m_obs'] = table.MaskedColumn(0, dtype=int, length=n)
    cat['n_g130m_obs'] = table.MaskedColumn(0, dtype=int, length=n)
    cat['lya_observed'] = table.MaskedColumn(False, dtype=bool, length=n)
    cat['lya_transit_observed'] = table.MaskedColumn(False, dtype=bool, length=n)

    for i, target in enumerate(cat):
        name = target['pl_name']
        if np.ma.is_masked(name):
            name = 'TOI-{}'.format(target['toi'])
        targcoord = coord.SkyCoord(target['ra']*u.deg, target['dec']*u.deg)
        dist = obs_coords.separation(targcoord)
        matches = dist < 10*u.arcmin
        matchobs = obs[matches]

        # determine if Lya line has been observed
        spec, wave = matchobs['spec'], matchobs['wave']
        e140m_lya = (np.char.count(spec, 'E140M') > 0) & (wave < 2000)
        g140m_lya = (np.char.count(spec, 'G140M') > 0) & ((wave == 1222) | (wave == 1218))
        g130m_lya = ((np.char.count(spec, 'G130M') > 0)
                     & ((1215 > wave - 150) & (1215 < wave - 7))
                     | ((1215 > wave + 7) & (1215 < wave + 150)))

        lya_mask = e140m_lya | g140m_lya | g130m_lya
        lya_observed = np.sum(lya_mask) > 1
        cat['lya_observed'][i] = lya_observed
        cat['n_e140m_obs'][i] = np.sum(e140m_lya)
        cat['n_g140m_obs'][i] = np.sum(g140m_lya)
        cat['n_g130m_obs'][i] = np.sum(g130m_lya)
        pids = np.unique(matchobs[lya_mask]['prop'])
        obsabs = '\n'.join([abdic[pid] for pid in pids])
        cat['hst abstracts'][i] = obsabs

        if lya_observed:
            # I was about to just do the math to see if any observations align with transit, but this has two problems:
            # (1) the duplication database doesn't have data on when the observation occurred, so I'd have to pull in another databse
            # (2) won't work for planned observations not yet executed

            print('')
            print(name)
            print('------')
            print(obsabs)
            cols = 'config mode spec aper time prop dataset'.split()
            lyaobs = matchobs[lya_mask]
            lyaobs.sort(['prop', 'dataset'])
            lyaobs[cols].pprint(-1)
            print('\n\n\n')
            answer = input("Hit enter if lya obs are good, any letter if not.".format(name))
            if answer != '':
                cat['lya_observed'][i] = False
                cat['n_e140m_obs'][i] = 0
                cat['n_g140m_obs'][i] = 0
                cat['n_g130m_obs'][i] = 0
            answer = input("Mark transit as observed for {} (y/n)?".format(name))
            cat['lya_transit_observed'][i] = answer == 'y'
        else:
            print('No lya observations for {}.'.format(name))
    # endregion

    cols = ['pl_name', 'toi', 'ra', 'dec']
    lya_observed = cat[cols][cat['lya_observed']]
    transit_observed = cat[cols][cat['lya_transit_observed']]
    return lya_observed, transit_observed


def filter_TESS_candidates(combo):
    x = combo['tfopwg_disp'].filled('PC')
    return ((x == 'PC') | (x == 'APC')) & ~combo['toi'].mask


def count_stars(cat):
    Nra = len(np.unique(cat['ra']))
    Ndec = len(np.unique(cat['dec']))
    assert Nra == Ndec
    return Nra


# from this point on the rest is done by hand

#region Lya transit analysis
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


def ism_checker(systems,
                n_H0=0.06/u.cm**3,
                n_H1=0.01/u.cm**3,
                rows=3, cols=4, figsize=[11, 5.5]):
    """systems should be some portion of the "combo" table referenced in all
    the functions above that has systems I want to investigate more closely"""
    Tism = 10000 * u.K
    wgrid = np.arange(1210., 1220., 0.005) * u.AA
    exp = 14

    etc_flux = 1e-14
    etc_src = 105*2/10000.
    etc_bg = 10.5*2/10000.

    panels = rows*cols
    for i in range(0, len(systems), panels):
        fig, axs = plt.subplots(rows, cols, figsize=figsize, dpi=100)
        fig.subplots_adjust(0.07, 0.07, 0.94, 0.9, hspace=0.05, wspace=0.2)
        bax = utils.common_axes(fig)
        bax2 = bax.twinx()
        [bax2.spines[s].set_visible(False) for s in
         ['top', 'bottom', 'left', 'right']]
        bax2.tick_params(labelleft=False, labelbottom=False, labelright=False,
                         left='off', bottom='off', right='off')
        bax2.set_zorder(-10)
        bax.tick_params(left='off')
        bax.set_xlabel('Velocity Relative to Planet [km/s]', labelpad=16)
        bax.set_ylabel('Flux at Earth [$10^{{-{}}}\ \mathrm{{erg\ s^{{-1}}\ '
                       'cm^{{-2}}\ \AA^{{-1}}}}$]'.format(exp), labelpad=25)
        bax2.set_ylabel('ISM Transmittance', labelpad=25)

        for i, (star, ax) in enumerate(zip(systems[i:i+panels], axs.ravel())):
            name = star['pl_name']
            if not name:
                name = star['toi']
            if not name:
                name = star['epic_candname']
            ax2 = ax.twinx()
            ax2.tick_params(labelright=False)
            sim_file = glob.glob('simulation results/*{}*.csv'.format(name))
            if len(sim_file) == 1:
                data = np.loadtxt(sim_file[0], delimiter=',', skiprows=1)
                v_trans, _, d_trans, _ = data.T
                v_trans = v_trans/1e5
                max_depth = np.max(d_trans)
            else:
                v_trans, d_trans = np.loadtxt(paths.reffolder / 'generic profile.csv', delimiter=',').T
                # scale depth to match extended hill sphere transit
                max_depth = depth_extended_hill_transit(star)
                d_trans = max_depth / np.max(d_trans) * d_trans
            for n_H in (n_H1, n_H0):
                Nh = n_H*star['sy_dist'] * u.pc
                rv_star = star['st_radv'] * u.km/u.s
                rv_ism = ism_velocity(star['ra'], star['dec'])
                lya = star['Flya_at_earth'] * u.Unit('erg s-1 cm-2')
                v = (wgrid/wlab_H - 1)*3e5 - rv_star.value

                intrinsic = reversed_lya_profile(wgrid, rv_star, lya)*10**exp
                transmitted = transmission(wgrid, rv_ism, Nh, Tism)
                observed = intrinsic * transmitted
                transit_depth = np.interp(v, v_trans, d_trans)
                transit = observed * (1 - transit_depth)

                intr, = ax.plot(v, intrinsic, color='C0')
                obs, = ax.plot(v, observed, color='C1')
                _ = ax.fill_between(v, transit, observed, color='C1',
                                    alpha=0.4)
                tt, = ax.plot(v,  transit, color='C1', ls=':')
                tr, = ax2.plot(v, transmitted, color='k', alpha=0.5, lw=1)

                if n_H == n_H0:
                    fmx = np.max(intrinsic).value/2
                    src_cps = fmx*10**-exp/etc_flux * etc_src
                    src_cnts = src_cps* 3 * 2700.
                    bg_cnts = etc_bg * 3 * 2700.
                    err = np.sqrt(src_cnts + bg_cnts)/src_cnts*fmx
                    ax.errorbar(175, fmx, err, fmt='.', color='k', capsize=1)

            if (i/2 < rows - 1) and (i < panels - 1):
                ax.tick_params(labelbottom=False)

            ax.set_xlim(-200, 200)
            _,y1 = ax.get_ylim()
            ax.set_ylim(None, y1*1.15)
            name = ('{}  {:.1f} km/s  {:.0f} pc  {:.2f}%'
                    '\n{:.0f} K'
                    ''.format(name, rv_star.value, star['sy_dist'], max_depth*100, star['pl_eqt']))
            ax.annotate(name, xy=(0.5, 0.96), xycoords='axes fraction',
                        bbox=dict(fc='w', lw=0, alpha=0.3, pad=2),
                        ha='center', va='top', fontsize='small')
            ax2.set_ylim(-0.05, 1.15)
            ax2.set_yticks((0, 0.5, 1))

        fig.legend((intr, obs, tt, tr), ('Intrinsic Ly$\\alpha$',
                                     'Observed Ly$\\alpha$',
                                         'In-Transit Ly$\\alpha$',
                                         'ISM Transmittance'),
                   loc='upper center', ncol=4)


ism_columns = table.Table.read(paths.reffolder / 'redfield_N_H_columns.ecsv')
n_H = 10 ** ism_columns['N_H'] * u.cm ** -2 / (ism_columns['d'] * u.pc)
n_H = n_H.to('cm-3')
ism_columns['n_H'] = n_H

def transit_snr(planet, expt_out=3500, expt_in=6000, add_jitter=False, optimistic=False, pessimistic=False):
    if (optimistic and add_jitter) or (pessimistic and add_jitter):
        raise ValueError
    wgrid = np.arange(1210., 1220., 0.005) * u.AA

    def jitter(value, err):
        return np.random.randn(1)*err + value
    def jitter_dex(value, errdex):
        return value*10**(np.random.randn(1) * errdex)

    if add_jitter:
        n_H = np.random.choice(ism_columns['n_H'])
        lya_factor = jitter_dex(1, 0.2)
    else:
        n_H = np.median(ism_columns['n_H'])
        lya_factor = 1
    if optimistic:
        n_H = np.percentile(ism_columns['n_H'], 16)
        lya_factor = 10**0.2
    if pessimistic:
        n_H = np.percentile(ism_columns['n_H'], 95)
        lya_factor = 10**-0.4
    n_H /= u.cm**3

    v_integrate = [-150, 100]

    rv_star = planet['st_radv'] * u.km / u.s
    v_trans, d_trans_temp = np.loadtxt(paths.reffolder / 'generic profile.csv', delimiter=',').T
    # broaden the absorption profile to better match GJ 436 transit
    v_trans *= 2
    d_trans_norm = d_trans_temp / np.max(d_trans_temp)
    v = (wgrid / lyasim.wlab_H - 1) * 3e5 - rv_star.value
    v_etc_PF = lyasim.v_etc - rv_star.value

    # RVing the planet parameters, if needed
    if add_jitter:
        planet = copy(planet)
        def jittercol(col, default_dex):
            if np.ma.is_masked(planet[col+'err1']):
                planet[col] = jitter_dex(planet[col], default_dex)
            else:
                planet[col] = jitter(planet[col], planet[col+'err1'])
        rdex = 0.075 if planet['pl_rade'] > 1.23 else 0.045
        jittercol('pl_bmasse', rdex)
        jittercol('st_mass', 0.04)
        jittercol('pl_orbsmax', 0.04)
        jittercol('st_rad', 0.04)

    max_depth = lyasim.depth_extended_hill_transit(planet)
    d_trans = max_depth * d_trans_norm
    observed = lyasim.lya_at_earth_auto(planet, n_H, lya_factor)
    v_blue = v[np.argmax(observed[v < 0])]
    v_trans_shifted = v_trans + v_blue
    transit_depth = np.interp(v, v_trans_shifted, d_trans)
    intransit = observed * (1 - transit_depth)
    fo, eo = lyasim.sim_g140m_obs(observed, expt_out)
    fi, ei = lyasim.sim_g140m_obs(intransit, expt_in)
    d = fo - fi
    e = np.sqrt(eo ** 2 + ei ** 2)

    integrate = (v_etc_PF > v_integrate[0]) & (v_etc_PF < v_integrate[1])
    D = np.sum(d[integrate])
    E = np.sqrt(np.sum(e[integrate]**2))
    return D/E
# endregion
