import re

column_map = {
    'st_plx': 'sy_plx',
    'st_dist': 'sy_dist',
    'st_pmra': 'sy_pmra',
    'st_pmdec': 'sy_pmdec',
    'st_tmag': 'sy_tmag',
    'pl_trandurh': 'pl_trandur',
    'st_kep': 'sy_kepmag',
    'st_bj': 'sy_bmag',
    'st_vj': 'sy_vmag',
    'st_us': 'sy_umag',
    'st_gs': 'sy_gmag',
    'st_rs': 'sy_rmag',
    'st_is': 'sy_imag',
    'st_zs': 'sy_zmag',
    'st_j2': 'sy_jmag',
    'st_h2': 'sy_hmag',
    'st_k2': 'sy_kmag',
    'st_wise1': 'sy_w1mag',
    'st_wise2': 'sy_w2mag',
    'st_wise3': 'sy_w3mag',
    'st_wise4': 'sy_w4mag',
    'st_vsini': 'st_vsin',
    'st_metfe': 'st_met',
    'ra_str': 'rastr',
    'dec_str': 'decstr',
    'pl_pnum': 'sy_pnum',
    'tid': 'tic_id'
}


useful_columns = [
 'rowid',
 'pl_name',
 'hostname',
 'pl_letter',
 'hd_name',
 'hip_name',
 'tic_id',
 'gaia_id',
 'cb_flag',
 'discoverymethod',
 'disc_year',
 'disc_refname',
 'disc_pubdate',
 'rv_flag',
 'tran_flag',
 'pl_controv_flag',
 'pl_orbper',
 'pl_orbsmax',
 'pl_orbeccen',
 'pl_rade',
 'pl_bmasse',
 'pl_bmassprov',
 'pl_dens',
 'pl_insol',
 'pl_eqt',
 'pl_tranmid',
 'pl_tranmid_systemref',
 'pl_imppar',
 'pl_trandep',
 'pl_trandur',
 'pl_ratdor',
 'pl_ratror',
 'pl_rvamp',
 'st_spectype',
 'st_teff',
 'st_rad',
 'st_mass',
 'st_met',
 'st_metratio',
 'st_lum',
 'st_logg',
 'st_age',
 'st_dens',
 'st_vsin',
 'st_rotp',
 'st_radv',
 'rastr',
 'ra',
 'decstr',
 'dec',
 'sy_dist',
 'sy_plx',
 'sy_pmra',
 'sy_pmdec',
 'sy_bmag',
 'sy_vmag',
 'sy_jmag',
 'sy_hmag',
 'sy_kmag',
 'sy_umag',
 'sy_gmag',
 'sy_rmag',
 'sy_imag',
 'sy_zmag',
 'sy_w1mag',
 'sy_w2mag',
 'sy_w3mag',
 'sy_w4mag',
 'sy_gaiamag',
 'sy_icmag',
 'sy_tmag',
 'sy_kepmag',
 'pl_nnotes',
 'epic_name',
 'tm_name',
 'epic_candname',
 'k2c_refdisp',
 'k2c_disp',
 'k2c_note',
 'k2_campaign',
 'k2c_recentflag',
 'pl_fppprob',
 'toi',
 'toipfx',
 'tid',
 'ctoi_alias',
 'tfopwg_disp',
 'toi_created',
 'rowupdate',
 'id',
 'disposition'
]


def homogenize_columns(cat):
    """Rename columns of K2 or TESS exoplanet candidate catalogs to match the
    composite system parameters catalog. """

    # deal with possibility of only symerr or err
    for name in cat.colnames:
        if name in column_map:
            newname = column_map[name]
            cat.rename_column(name, newname)
            has_dual_err = name + 'err1' in cat.colnames
            if has_dual_err:
                cat.rename_column(name + 'err1', newname + 'err1')
                cat.rename_column(name + 'err2', newname + 'err2')
            if name + 'err' in cat.colnames:
                if not has_dual_err:
                    cat[newname + 'err1'] = cat[name + 'err']
                    cat[newname + 'err2'] = cat[name + 'err']
                cat.remove_column(name + 'err')
            if name + 'symerr' in cat.colnames:
                if not has_dual_err:
                    cat[newname + 'err1'] = cat[name + 'symerr']
                    cat[newname + 'err2'] = cat[name + 'symerr']
                cat.remove_column(name + 'symerr')


def basename(name):
    base = re.sub('(sym)?err[12]?', '', name)
    base = base.replace('reflink', '')
    base = base.replace('lim', '')
    base = base.replace('src', '')
    return base


def remove_dispensible_columns(cat):
    for name in cat.colnames:
        base = basename(name)
        if base not in useful_columns:
            cat.remove_column(name)