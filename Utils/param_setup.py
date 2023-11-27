from astropy.io import ascii
import random
import glob
from astropy.table import Table, join


def SELFIE_extractor() -> Table:

    """
    Extracts all relevant parameters from set of simulation results.

    :return: Table of parameters.
    :rtype: Table
    """

    # import all of our tables
    print('loading gal_2pt0')
    gal_2pt0 = Table.read("/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/"
                          "fits_storage_2/2pt0_wfd_galdat.csv",
                          format='csv', delimiter=',')
    print('gal_2pt0 loaded successfully')
    print(len(gal_2pt0))
    print('loading phase_data')
    phase_data = Table.read("/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/"
                            "Fits_storage_2/SNe_phasedata1.csv",
                            format='csv', delimiter=',')
    print('phase data loadded successfully')
    print('loading SELFIEsim_SNe')
    SELFIE_SNe = Table.read("/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/"
                            "Fits_storage_2/SELFIE172_S1001.csv",
                            format='csv', delimiter=',')
    print('SELFIEsim_SNe loaded successfully')
    print('loading cat_export_SNe')
    cat_export_SNe = Table.read("/Users/andrew/Desktop/Python_Stuff/"
                                "SN_and_Galaxy/Fits_storage_2/"
                                "catex45_S1001.csv",
                                format='csv', delimiter=',')
    print('cat_export_SNe loaded successfully')

    print('loading texp data')
    tiles = Table.read("/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/"
                       "fits_storage/SELFIE172_tiles.fits")
    fibres = Table.read("/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/"
                        "fits_storage/SELFIE172_fibres.fits")

    extracted_date_info = join(fibres['targ_id', 'tile_id'],
                               tiles['tile_id', 'texp', 'jd_obs'],
                               keys='tile_id', join_type='inner')

    # must start by linking to cat_export_SNe as this has a name column.
    # Cut down to only be the useful columns at the rows where
    # the name value is lower than the max present in the 2.0 galaxy data

    new_table = join(cat_export_SNe['ra', 'dec', 'mag', 'redshift_estimate',
                                    'targ_id', 'u_obj_id', 'name'],
                     gal_2pt0['redshift_final', 'hostgal_mag_r', 'ra', 'dec',
                              'sim_model_name', 'peakmjd', 'sim_type_name',
                              'hostgal_snsep', 'name', 'hostgal_ddlr'],
                     keys='name', join_type='inner')
    print('galaxy 2pt0 and catalog export linked successfully')

    # mask the SELFIE data so only fobs > whatever are accepted
    mask = (SELFIE_SNe['fobs'] != 0)
    SELFIE_SNe = SELFIE_SNe[mask]

    # join the table to the SEFLIE data
    new_table = join(new_table,
                     SELFIE_SNe['ra', 'dec', 'u_obj_id', 'fobs',
                                'jd_obs_first', 'jd_obs_last',
                                'jd_last_observed', 'texp_d',
                                'texp_g'],
                     keys='u_obj_id', join_type='inner')
    print('SELFIE sim linked successfully')

    # now join to the phase data
    new_table = join(new_table,
                     phase_data['name', 'RA', 'DEC', 'mag',
                                'redshift_estimate', 'TEMPLATE'],
                     keys=['mag', 'redshift_estimate', 'name'],
                     join_type='inner')
    print('phase data linked successfully')

    # add the texp (this is from the product file of date extractor,
    # so may need to adjust this to get the raw process)
    new_table = join(new_table,
                     extracted_date_info['texp', 'targ_id', 'jd_obs'],
                     keys=['targ_id'])
    print('date data linked successfully')

    # this seems to find some duplicate rows (every row except TEMPLATE)
    # so remove them here

    df_state = Table.to_pandas(new_table)

    df_nodupe = df_state.drop_duplicates(
        subset=['ra_1', 'dec_1', 'mag', 'redshift_estimate', 'targ_id',
                'u_obj_id', 'name', 'redshift_final', 'hostgal_mag_r',
                'ra_2', 'dec_2', 'sim_model_name', 'peakmjd',
                'sim_type_name', 'hostgal_snsep', 'ra',
                'dec', 'fobs', 'jd_obs_first', 'jd_obs_last',
                'jd_last_observed', 'RA', 'DEC'], keep='last')

    to_save = Table.from_pandas(df_nodupe)
    print(len(df_state), len(df_nodupe))

    print(df_nodupe.columns)

    # still need to get the actual phase value
    phase_val = []
    for c in range(len(to_save)):
        start = str(to_save[c][26]).find('phase') + 5
        end = str(to_save[c][26]).find('_red')
        phase_val.append(float(to_save[c][26][start:end]))

    # add the extra row to the table and then save and return it
    to_save['phase_val'] = phase_val

    # generate exposure times
    texp_obj = []
    for i in range(len(to_save)):
        texp_obj.append(to_save['texp_d'][i] * to_save['fobs'][i])

    to_save['texp_obj'] = texp_obj

    print('first line of table being saved = ', to_save[0])

    ascii.write(to_save,
                "/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/Results/"
                "text_dump_SNANAtemps/SELFIE172_SNANA_tests_WFD.txt",
                format='csv', delimiter=',', overwrite=True)

    return to_save


def adj_setup(data: Table) -> list:

    """
    Takes table of parameters extracted from simulations and
    adds then to lists for use in generating spectra.

    :param table data: Table of parameters.
    :return: list of values to loop into comb_maker().
    :rtype: list
    """

    Smags = []
    Gmags = []
    redshift = []
    templates = []
    SNe_types = []
    phase = []
    ddlr = []
    snsep = []
    texp_visit = []
    texp_obj = []

    real_data = data

    for it in range(len(real_data)):
        Smags.append(real_data['mag'][it])
        Gmags.append(real_data['hostgal_mag_r'][it])
        redshift.append(real_data['redshift_estimate'][it])
        templates.append(real_data['TEMPLATE'][it])
        SNe_types.append(real_data['sim_type_name'][it])
        phase.append(-99)
        ddlr.append(real_data['hostgal_ddlr'][it])
        snsep.append(real_data['hostgal_snsep'][it])
        texp_visit.append(38.712673)
        texp_obj.append(26.323857223591407)

    gal_array = glob.glob(
        '/Users/Andrew/Desktop/Python_Stuff/SN_and_Galaxy/galaxy_templates/*'
        )
    SNe_temp_array = glob.glob(
        '/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/SNANA_temps/*'
        )

    SNe_temp_names = []
    for t in range(len(SNe_temp_array)):
        SNe_name_begins = SNe_temp_array[t].find('SNANA_temps/') + 12
        SNe_temp_names.append(SNe_temp_array[t][SNe_name_begins:])

    galaxies = []
    supernovae = []
    for ti in range(len(Smags)):

        galaxies.append(random.choice(gal_array))
        supernovae.append(SNe_temp_array[SNe_temp_names.index(templates[ti])])

    return [redshift, Smags, Gmags, phase, templates, SNe_types, galaxies,
            supernovae, ddlr, snsep, texp_obj, texp_visit]
