from pathlib import Path

from astropy import table
from astropy import units as u
from astropy import coordinates as coord
import numpy as np
from tqdm import tqdm
from scipy import interpolate

import finder_selecter as fs
import galex_motion
import lyasim
import paths
import utilities as utils


#%% *NEEDED DOWNSTREAM* load in the stage 1 target list and filter to unique stars

targets = table.Table.read(paths.stage1_targets)
targets = fs.unique_stars(targets)
n = len(targets)

# delete misleading, vestigial transit SNR columns (numbers are very wrong)
targets.remove_columns(['pl_lya_tnst_snr_nominal', 'pl_lya_tnst_snr_optimistic'])


#%% *NEEDED DOWNSTREAM* swap out bad targets for next best from stage 1 menu
remove = ('WD 1856+534', 'WASP-189')
# A star, white dwarf

prelim = table.Table.read(paths.preliminary_targets)

longperiod = prelim['pl_orbper'] > 90
FPs = np.zeros(len(prelim), bool)
toi_false_positives = np.loadtxt(paths.toi_false_positives)
for toi in toi_false_positives:
    match = prelim['toi'].filled('') == str(toi)
    if not np.any(match):
        print('No match found for {}'.format(toi))
    FPs = FPs | match
noobs = ~prelim['lya_observed'].data
cool = prelim['st_teff'] < 6500

# compute crude SNRs to compare with ethan's estimates
prelim['lya_transit_snr_crude'] = 0.0
prelim['lya_transit_snr_crude_optimistic'] = 0.0
for i, planet in enumerate(prelim):
    prelim['lya_transit_snr_crude'][i] = fs.transit_snr(planet)
    prelim['lya_transit_snr_crude_optimistic'][i] = fs.transit_snr(planet, optimistic=True)

# flag targets that were previously cut so I can compare with new cuts
special_cases = ~prelim['flag_snr_okay'].data
bad_initial_SNR = ((prelim['lya_transit_snr_crude_optimistic'] < 2)
                   | ((prelim['lya_transit_snr_crude'] < 2) & ~special_cases))
SNRratio = prelim['lya_transit_snr_sim'] / prelim['lya_transit_snr_crude']
bad_SNR_estimate = ((SNRratio < 0.1)
                    & (SNRratio > 1e-10))
usable = ~longperiod & ~FPs & ~bad_SNR_estimate & ~bad_initial_SNR
oldcuts = ~(usable & noobs)
prelim['cut_from_proposal'] = oldcuts

# now make a new cut that differs from previous cut in order to provide more options
notransobs = ~prelim['lya_transit_observed'] # GJ 1132 leaks in otherwise
detectable = (prelim['lya_transit_snr_sim_optimistic'] > 3) | (prelim['lya_transit_snr_crude_optimistic'] > 3)
usable = ~longperiod & ~FPs & detectable & cool
keep = usable & noobs
prelim = prelim[keep]

# this is so you can have a look at differences
# prelim['crude/sim'] = prelim['lya_transit_snr_crude_optimistic']/prelim['lya_transit_snr_sim_optimistic']
# prelim['crude/sim'][prelim['lya_transit_snr_sim_optimistic'] < 0.1] = 0
# prelim.sort('crude/sim', reverse=True)
# prelim['Enom'] = prelim['lya_transit_snr_sim']
# prelim['Eopt'] = prelim['lya_transit_snr_sim_optimistic']
# prelim['Cnom'] = prelim['lya_transit_snr_crude']
# prelim['Copt'] = prelim['lya_transit_snr_crude_optimistic']
# prelim['Cnom'].format='.1f'
# prelim['Copt'].format='.1f'
# prelim['Enom'].format='.1f'
# prelim['Eopt'].format='.1f'
# prelim['crude/sim'].format='.1f'
# showcols = 'hostname st_spectype pl_rade Cnom Copt Enom Eopt cut_from_proposal crude/sim'.split()
# prelim[prelim['cut_from_proposal']][:10][showcols].pprint(-1,-1)

# make an even shorter list of stars to choose from
prelim['crude/sim'] = prelim['lya_transit_snr_crude_optimistic']/prelim['lya_transit_snr_sim_optimistic']
prelim['crude/sim'][prelim['lya_transit_snr_sim_optimistic'] < 0.1] = 0
keep = ((prelim['lya_transit_snr_sim_optimistic'] > 3)
        | ((prelim['crude/sim'] > 5) & (prelim['lya_transit_snr_sim_optimistic'] > 1)))
keep = keep & prelim['flag_abovegap']
shortlist = prelim[keep]

shortlist.sort('lya_transit_snr_sim_optimistic', reverse=True) # I know this is a duplicate line but it is critical to the loop so want to be sure it gets run even if I delete eearlier code
i = 0
while len(targets) < n + len(remove):
    row = shortlist[i]
    ra_match = np.isclose(targets['ra'], row['ra'])
    if not np.any(ra_match):
        targets.add_row(None)
        for col in targets.colnames:
            if col in shortlist.colnames:
                if targets[col].dtype != shortlist[col].dtype: # fixes issues with differing string lengths
                    targets[col] = targets[col].astype(shortlist[col].dtype)
                targets[col][-1] = shortlist[col][i]
    i = i + 1
for name in remove:
    mask = targets['hostname'] != name
    targets = targets[mask]


#%% *NEEDED DOWNSTREAM* fix the effective temperature of EPIC 205530323
# not sure why it is off. other effective temperatures in spot checks seem fine
# seems like maybe the effective temp changed in an update to the TIC bc other stars with only a TIC Teff
# seem to have changed some
i = targets['hostname'] == 'EPIC 205530323'
targets['st_teff'][i] = 3142.0


#%% *NEEDED DOWNSTREAM* query SIMBAD for info needed for APT
# need to replace ra, dec, and pms with results from simbad bc I don't trust exoplanet archive
# simbad recognizes TIC IDs and there is a TIC ID for every object
from astroquery.simbad import Simbad
simquery = Simbad()
fields = 'coo_err_maja plx pmra pmdec flux(V) flux_error(V) flux(G) flux_error(G) sptype sp_qual'.split()
simquery.add_votable_fields(*fields)
simbad = simquery.query_objects(targets['tic_id'])


#%% *NEEDED DOWNSTREAM* copy in spectral types where quality is better than C (this actually replaces nothing, but it does add src column)
s1_bad = targets['st_spectype'].mask | (targets['st_spectype'].filled('') == '')
simbad_good = (simbad['SP_QUAL'] <= 'C') & (simbad['SP_TYPE'] != '')
replace = s1_bad & simbad_good # this turns out to be zero :)
targets['st_spectype'][replace] = simbad['SP_TYPE'][replace]
fs.add_src_col(targets, 'st_spectype')
targets['st_spectypesrc'][~s1_bad] = 'NASA Exoplanet Archive'
targets['st_spectypesrc'][replace] = 'SIMBAD'



#%% *NEEDED DOWNSTREAM* get spectral types based on Teff for the remainder using mamajek table
mamajek = table.Table.read(paths.reffolder / 'mamajek_SpT_table.ecsv')
interper = interpolate.interp1d(mamajek['Teff'], np.arange(len(mamajek)), 'nearest')
bad = targets['st_spectype'].filled('') == ''
itype = interper(targets['st_teff'][bad])
mamSpT = mamajek['SpT'][itype.astype(int)]
targets['st_spectype'][bad] = mamSpT
targets['st_spectypesrc'][bad] = 'Mamajek table'


#%% *NEEDED DOWNSTREAM* update and fill in V mags for ETC calcs

# use latest V mag from simbad wherever possible
move = ~simbad['FLUX_V'].mask
targets['sy_vmag'][move] = simbad['FLUX_V'][move]
targets['sy_vmagerr1'][move] = targets['sy_vmagerr2'][move] = simbad['FLUX_ERROR_V'][move]


# copy over G mags from simbad
move = ~simbad['FLUX_G'].mask
targets['sy_gmag'][move] = simbad['FLUX_G'][move]

# use G mags to estimate V mag based on mamajek table where there isn't a V
assert not hasattr(targets['st_teff'], 'mask') # if this isn't true, need to deal with masked values
G_V = utils.safe_interp_table(targets['st_teff'], 'Teff', 'G-V', mamajek)
V_from_G = targets['sy_gmag'] - G_V
move = ~V_from_G.mask & targets['sy_vmag'].mask
targets['sy_vmag'][move] = V_from_G[move]

# this still leaves at least one star without a V mag. the only mag it has is a TESS mag.
# I'll estimate from the mamajek table based on Teff and distance.
Mv = utils.safe_interp_table(targets['st_teff'], 'Teff', 'Mv', mamajek)
d = targets['sy_dist'].quantity.to_value('pc')
V = Mv + 5*np.log10(d) - 5
move = targets['sy_vmag'].mask
targets['sy_vmag'][move] = V[move]

#%% save the Ms to the cycle 32 folder
Ms = np.char.startswith(targets['st_spectype'], 'M')
Ms = targets[Ms]
Ms.write(paths.output32 / 'M stars.ecsv')


#%% make APT-ready names and start APT table
names = targets['hostname'].copy()
mask = names.mask
names[mask] = np.char.add('TOI-', targets['toipfx'][mask]) # make sure aptcat has no indices or this will fail
names = np.char.replace(names, ' ', '')
aptcat = table.Table([names], names=['Target Name'], masked=True)


#%% APT sky coordinates
aptcat['RA'] = simbad['RA']
aptcat['RA Uncertainty'] = simbad['COO_ERR_MAJA']/1000 # this will be an upper limit, needs to be in arcsec
aptcat['DEC'] = simbad['DEC']
aptcat['DEC Uncertainty'] = simbad['COO_ERR_MAJA']/1000 # needs to be in arcsec
aptcat['Reference Frame'] = 'ICRS'
aptcat['Epoch'] = 2000.0

#%% APT proper motion and parallax
aptcat['RA PM'] = simbad['PMRA']
aptcat['RA PM units'] = str(simbad['PMRA'].unit).upper().replace(' ', '')
aptcat['DEC PM'] = simbad['PMDEC']
aptcat['DEC PM units'] = str(simbad['PMDEC'].unit).upper().replace(' ', '')
aptcat['Annual Parallax'] = simbad['PLX_VALUE']/1000 # needs to be in arcsec

#%% APT magnitudes
aptcat['V-Magnitude'] = targets['sy_vmag']
aptcat['Mag Uncertainty'] = targets['sy_vmagerr1']
Gfluxes = []
for i in range(n):
    if simbad['FLUX_G'].mask[i]:
        Gfluxes.append('')
    else:
        G = simbad['FLUX_G'][i]
        Gfluxes.append(f'G={G}')
aptcat['Other Fluxes'] = Gfluxes

#%% APT description
aptcat['Category'] = 'STAR'
descriptions = []
for i in range(n):
    SpT = targets['st_spectype'][i]
    SpT = SpT.replace(' ', '')
    assert 'I' not in SpT
    if '/' in SpT:
        SpT = SpT.split('/')[0]
    letter = SpT[0]
    if letter == 'F':
        num = 1 if len(SpT) == 1 else int(SpT[1])
        if num < 3:
            SpT = 'F0-F2'
        else:
            SpT = 'F3-F9'
    else:
        SpT = f'{letter} V-IV'
        desc = f'[{SpT}, Extra-solar Planetary System]'
    descriptions.append(desc)
aptcat['Description'] = descriptions
aptcat['Radial Velocity'] = targets['st_radv']


#%% APT save targets for APT input

# make sure you didn't skip any code blocks (as you have done once already :)
cols_from_blocks = 'Target Name,RA,RA PM,V-Magnitude,Category'.split(',')
for col in cols_from_blocks:
    assert col in aptcat.colnames

path = paths.output32 / 'targets_for_APT.csv'
aptcat.write(path, format='ascii.commented_header', overwrite=True, delimiter=',')



#%% get a median V mag to use in ACQ ETC for each Teff

pickles_Teffs = np.array([2951, 3548, 3801, 4188, 4886, 5188, 5584, 6039, 6531, 6776, 7211])
mids = (pickles_Teffs[1:] + pickles_Teffs[:-1])/2
bins = np.insert(mids, [0, len(mids)], [0, 10000])
for Ta, Tb in zip(bins[:-1], bins[1:]):
    inbin = (targets['st_teff'] >= Ta) & (targets['st_teff'] < Tb)
    Vmin = np.min(targets['sy_vmag'][inbin])
    Vmed = np.median(targets['sy_vmag'][inbin])
    Vmax = np.max(targets['sy_vmag'][inbin])
    print(f'Teff = {Ta}—{Tb} min-med-max V {Vmin:.1f}–{Vmed:.1f}–{Vmax:.1f}')


#%% *NEEDED DOWNSTREAM* clear out any old galex mags
for col in targets.colnames:
    if 'uvmag' in col:
        targets[col].mask = True


#%% *NEEDED DOWNSTREAM* pull galex mags. redo this a few times to ensure any failed requests get another chance.

# gotta use the simbad coordinates. archive coordinates are at a weird equinox

def add_mag(cps, cpserr, band, i):
    if np.isnan(cps):
        pass
    elif cps < 0:
        result = galex_motion.galex_cps2mag(cpserr, band)
        if not np.isnan(result):
            targets[f'sy_{band}mag'][i] = result
            targets[f'sy_{band}maglim'][i] = -1
    else:
        result = galex_motion.galex_cps2mag(cps, band)
        if not np.isnan(result):
            targets[f'sy_{band}mag'][i] = result
            magerr = galex_motion.cps2magerr(cps, cpserr)
            targets[f'sy_{band}magerr1'][i] = magerr
            targets[f'sy_{band}magerr2'][i] = magerr
            targets[f'sy_{band}maglim'][i] = 0

redo_mask = targets['sy_nuvmag'].mask & targets['sy_fuvmag'].mask
redo_args, = np.nonzero(redo_mask)
for i in tqdm(redo_args):
    try:
        ra = simbad['RA'][i]
        dec = simbad['DEC'][i]
        coo = coord.SkyCoord(ra, dec, unit=('hourangle', 'deg'))
        ra, dec = coo.ra.value, coo.dec.value
        if simbad['PMRA'].mask[i]:
            pm_ra = 0
            pm_dec = 0
        else:
            pm_ra = simbad['PMRA'][i]
            pm_dec = simbad['PMDEC'][i]

        result = galex_motion.extract_and_coadd(ra, dec, pm_ra, pm_dec, match_radius=16. / 3600, query_timeout=60.)
        (nuv, nerr), (fuv, ferr) = result
        add_mag(nuv, nerr, 'nuv', i)
        add_mag(fuv, ferr, 'fuv', i)
    finally:
        continue

assert not np.any(np.isnan(targets['sy_fuvmag']))

#%% *NEEDED DOWNSTREAM* predict missing galex mags

# using a catalog of galex mags vs B-V color to do estimate a mag for a standard field star and then adjusting
# based on rotation rate or assuming saturated for stars without a known rotation rate or age

catsup = table.Table.read(paths.reffolder / 'catsup.vot')
usable = ~catsup['B-V'].mask & ~catsup['FUVmag'].mask
catsup = catsup[usable]
catsup['FUV-V'] = catsup['FUVmag'] - catsup['Vmag']
outlier = (catsup['FUV-V'] > 10) & (catsup['B-V'] < 0.25)
catsup = catsup[~outlier]
BV, FV = catsup['B-V'], catsup['FUV-V']
BVbreak = 0.8

# first fit to determine the "median"
lo = BV < BVbreak
pmlo = np.polyfit(BV[lo], FV[lo], 2)
pmhi = np.polyfit(BV[~lo], FV[~lo], 1)
def medfunc(BVvec):
    result = np.zeros_like(BVvec)
    lo = BVvec < BVbreak
    result[lo] = np.polyval(pmlo, BVvec[lo])
    result[~lo] = np.polyval(pmhi, BVvec[~lo])
    return result

# then fit to the upper half
upper = FV > medfunc(BV)
plo = np.polyfit(BV[upper & lo], FV[upper & lo], 2)
phi = np.polyfit(BV[upper & ~lo], FV[upper & ~lo], 1)
def upperfunc(BVvec):
    result = np.zeros_like(BVvec)
    lo = BVvec < BVbreak
    result[lo] = np.polyval(plo, BVvec[lo])
    result[~lo] = np.polyval(phi, BVvec[~lo])
    return result

BV = utils.safe_interp_table(targets['st_teff'], 'Teff', 'B-V', mamajek)
FV = upperfunc(BV)
FUVmag_field = FV + targets['sy_vmag']
saturated = targets['st_Ro'].filled(0) < 0.1 # assume saturated as a worst case
old = targets['st_age'].filled(0) > 1 # but not for stars we know are old
saturated[old] = False
Ro_sat = 0.1
Ro_nom = 0.7
index = 1.8
fluxfac = (targets['st_Ro'].filled(Ro_nom)/Ro_nom)**-index # note filling with Ro_nom here
fluxfac[saturated] = (Ro_sat/Ro_nom)**-index # then replace for stars we know or assume are saturatioed
mag_offset = -2.5*np.log10(fluxfac)
FUVmag = FUVmag_field - mag_offset

lims = targets['sy_fuvmaglim'].filled(0) == -1
replace = lims & (FUVmag.filled(0) > targets['sy_fuvmag'].filled(999))
targets['sy_fuvmag'][replace] = FUVmag[replace]

mask = targets['sy_fuvmag'].mask
targets['sy_fuvmag'][mask] = FUVmag[mask]


#%% add lya to the reference spectra (and fill in O I and negatives with zeros) and save

ref_spec_folder = Path(paths.refspecfolder / 'lya_masked')
etc_spec_folder = Path(paths.refspecfolder / 'lya_added_intrinsic')
files = list(ref_spec_folder.glob('*.spec'))
for path in files:
    spec = table.Table.read(path, format='ascii.ecsv')
    path_pieces = path.name.split('-')
    Teff = float(path_pieces[1])
    Prot = 1 if Teff < 2700 else 50 # this just puts the star in the right activity category for the linsky relationship
    Flya_1au, = fs.Lya_from_Teff_linsky13(np.array([Teff]), np.array([Prot]))
    d = spec.meta['distance']
    Flya_at_dist = Flya_1au * (u.AU/(d*u.pc))**2
    Flya_at_dist = Flya_at_dist.to_value('')
    lyamask = np.isnan(spec['f']) & (spec['w'] < 1250)
    wlya = spec['w'][lyamask]
    ylya = lyasim.reversed_lya_profile(wlya*u.AA, 0*u.km/u.s, Flya_at_dist*u.Unit('erg s-1 cm-2'))
    spec['f'][lyamask] = ylya.to_value('erg s-1 cm-2 AA-1')
    fillmask = np.isnan(spec['f']) | (spec['f'] < 0)
    spec['f'][fillmask] = 0

    savepath = etc_spec_folder / path.name.replace('.spec', '.dat')
    data = np.array((spec['w'].data, spec['f'].data))
    np.savetxt(savepath, data.T)


#%% get worst case fuv mag for each Teff bin of reference spectra to try in ETC

etc_spec_folder = Path('reference_uv_spectra/lya_added_intrinsic')
files = list(etc_spec_folder.glob('*.dat'))
Ts = []
for path in files:
    path_pieces = path.name.split('-')
    Teff = float(path_pieces[1])
    Ts.append(Teff)
Ts = np.sort(Ts)
Tmids = (Ts[:-1] + Ts[1:])/2
Tas = np.insert(Tmids, 0, 0)
Tbs = np.append(Tmids, 1e4)
for Ta, Tb in zip(Tas,Tbs):
    inbin = (targets['st_teff'] > Ta) & (targets['st_teff'] <= Tb)
    if sum(inbin) > 0:
        minfuv = np.nanmin(targets['sy_fuvmag'][inbin])
        print(f'{Ta}—{Tb} K min FUV {minfuv:.1f}')


#%% choose ACQ filter and add exposure times

etc_acq = table.Table.read(paths.etc32 / 'ACQ.csv')
targets['acq_filter'] = table.Column(length=n, dtype='object')
targets['acq_Texp_snr40'] = 0.0
targets['acq_Texp'] = 0.0
targets['acq_Tsat'] = 0.0

ietcs = list(range(len(etc_acq)))
findline = interpolate.interp1d(etc_acq['Teff'], ietcs, 'nearest')
T_hard_min = 0.1
LP = 'F28X50LP'
ND = 'F25ND3'
for i in range(n):
    Teff = targets['st_teff'][i]
    ietc = int(findline(Teff))
    Vtarg = targets['sy_vmag'][i]

    Texp_40 = {}
    Texp = {}
    Tsat = {}
    min_margin = {}
    for filter in [LP, ND]:
        Vetc = etc_acq[f'V_{filter}'][ietc]
        Texp_etc = etc_acq[f'Texp_{filter}'][ietc]
        Tsat_etc = etc_acq[f'Tsat_{filter}'][ietc]
        flux_ratio = 10**((Vetc - Vtarg)/2.5)
        Texp_40[filter] = Texp_etc / (flux_ratio)**2 # ^2 to get the same SNR as ETC
        Tsat[filter] = Tsat_etc / flux_ratio

        # if the saturation time is less than the min exposure, it's not going to work
        if Tsat[filter] < T_hard_min:
            Texp[filter] = np.nan
            min_margin[filter] = 0
            continue

        # try to split the difference between saturation and SNR 40 with a geometric mean
        Texp_desired = np.sqrt(Texp_40[filter]*Tsat[filter])
        if Texp_desired < T_hard_min:
            Texp_desired = T_hard_min
        Texp[filter] = Texp_desired
        bright_margin = Tsat[filter] / Texp_desired
        faint_margin = Texp_desired / Texp_40[filter]
        min_margin[filter] = min(faint_margin, bright_margin)

    # pick the filter that gives better margins
    # if both margins are > 10, pick the shorter exposure
    if (min_margin[LP] > 10) and (min_margin[ND] > 10):
        filter = LP
    elif min_margin[LP] > min_margin[ND]:
        filter = LP
    else:
        filter = ND
    targets['acq_filter'][i] = filter
    targets['acq_Texp_snr40'][i] = Texp_40[filter]
    targets['acq_Texp'][i] = Texp[filter]
    targets['acq_Tsat'][i] = Tsat[filter]


#%% estimate buffer times for G140L and G140M exposures

for grating in ['G140M','G140L']:
    etc = table.Table.read(paths.etc32 / f'{grating}.csv')
    ietcs = list(range(len(etc)))
    findline = interpolate.interp1d(etc['Teff'], ietcs, 'nearest', bounds_error=False, fill_value='extrapolate')
    Tbuf_col = f'Tbuffer_{grating}'
    targets[Tbuf_col] = 0.0

    for i in range(n):
        Teff = targets['st_teff'][i]
        ietc = int(findline(Teff))
        FUVtarg = targets['sy_fuvmag'][i]
        FUVetc = etc['FUVmag'][ietc]
        flux_ratio = 10 ** ((FUVetc - FUVtarg) / 2.5)
        Tbuffer_etc = etc['buffer_time'][ietc]
        Tbuffer = Tbuffer_etc / flux_ratio
        targets[Tbuf_col][i] = Tbuffer * 4/5 # 4/5 is based on recommendation on ETC



#%% print info for input to APT

# for pretty printing
targets['acq_Texp_snr40'].format = '.3g'
targets['acq_Texp'].format = '.3g'
targets['acq_Tsat'].format = '.3g'
targets['Tbuffer_G140M'].format = '.0f'
targets['Tbuffer_G140L'].format = '.0f'
targets['sy_vmag'].format = '.1f'

# cols = 'hostname toipfx st_teff sy_vmag acq_filter acq_Texp_snr40 acq_Texp acq_Tsat Tbuffer_G140M Tbuffer_G140L'.split()
cols = 'hostname toipfx acq_filter acq_Texp_snr40 acq_Texp acq_Tsat Tbuffer_G140M Tbuffer_G140L'.split()
targets[cols].pprint(-1,-1)

#%% save the final target catalog

targets.write(paths.output32 / 'cycle32_processed_targets.ecsv')