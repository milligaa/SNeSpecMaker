"""This module contains all of the utility functions that are
called when running the spec_maker module."""

from astropy.io import ascii
import random
import glob
from astropy.table import Table, join
import pandas as pd
from astropy.io import fits
import warnings

__docformat__ = 'reStructuredText'

def SELFIE_extractor(sim_file_dict: dict,
                     table_save_path: str) -> Table:

    """
    Extracts all relevant parameters from set of simulation results.

    :param sim_file_dict: dictionary containing paths to all sim
        files which can be produced by the assign_sim_files function.
    :param table_save_path: save table of combined
        data from all sims.
    :type sim_file_dict: dict
    :type table_save_path: str
    :return: Table of parameters.
    :rtype: Table.
    """

    # import all of our tables
    print('loading gal_2pt0')
    gal_2pt0 = Table.read(sim_file_dict["galpop"],
                          format='csv', delimiter=',')
    print('gal_2pt0 loaded successfully')
    print(len(gal_2pt0))
    print('loading phase_data')
    phase_data = Table.read(sim_file_dict["phase"],
                            format='csv', delimiter=',')
    print('phase data loadded successfully')
    print('loading SELFIEsim_SNe')
    SELFIE_SNe = Table.read(sim_file_dict["selfie"],
                            format='csv', delimiter=',')
    print('SELFIEsim_SNe loaded successfully')
    print('loading cat_export_SNe')
    cat_export_SNe = Table.read(sim_file_dict["catex"],
                                format='csv', delimiter=',')
    print('cat_export_SNe loaded successfully')

    print('loading texp data')
    tiles = Table.read(sim_file_dict["tiles"])
    fibres = Table.read(sim_file_dict["fibres"])

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
                (table_save_path+"SELFIE172_SNANA_tests_WFD.txt"),
                format='csv', delimiter=',', overwrite=True)

    return to_save


def adj_setup(data: Table) -> list:

    """
    Takes table of parameters extracted from simulations and
    adds then to lists for use in generating spectra.

    :param data: Table of parameters.
    :type data: Table
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
            supernovae, ddlr, snsep]


def assign_sim_files(sim_data_path):
    """
    Finds all simulation data files and determines which is which.

    :param sim_data_path: location of all sim files.
    :type sim_data_path: str
    """
    all_data_files = glob.glob((sim_data_path+'*'))
    print(all_data_files)

    selfie_file = ''
    catex_file = ''
    galpop_file = ''
    phasedata_file = ''
    fibres_file = ''
    tiles_file = ''

    selfie_cols = ['fobs', 'jd_obs_first', 'jd_obs_last',
                   'jd_last_observed', 'texp_d', 'texp_g']
    catex_cols = ['mag', 'targ_id', 'name']
    galpop_cols = ['sim_model_name', 'sim_type_name',
                   'hostgal_snsep', 'hostgal_ddlr']
    phase_cols = ['RA', 'DEC', 'TEMPLATE']
    tile_cols = ['texp', 'tile_id', 'jd_obs']
    fibre_cols = ['tile_id', 'targ_id']

    warnings.simplefilter("ignore")

    for file in all_data_files:
        print(file)
        try:
            hdul = fits.open(file, ignore_missing_simple=True)
            data_cols = list(hdul[1].data.columns)
            col_names = [col.name for col in data_cols]

            if all(x in col_names for x in tile_cols) is True:
                tiles_file = file

            elif all(x in col_names for x in fibre_cols) is True:
                fibres_file = file

            else:
                raise IndexError(("Correct fits columns not found "
                                  "check condition of "+file))

        except OSError:
            print(file, ' is not a fits file, retrying for csv')

            try:
                data_cols = pd.read_csv(file, index_col=None,
                                        nrows=0).columns.to_list()
                print(data_cols, '------------------------------------')

                if all(x in data_cols for x in selfie_cols) is True:
                    selfie_file = file

                elif all(x in data_cols for x in catex_cols) is True:
                    catex_file = file

                elif all(x in data_cols for x in galpop_cols) is True:
                    galpop_file = file

                elif all(x in data_cols for x in phase_cols) is True:
                    phasedata_file = file

                else:
                    raise IndexError((
                        'csv sim file columns not recognised for '+file))

            except UnicodeDecodeError or IndexError:
                raise TypeError((file+' is not a valid simulation file'))

    file_dict = {
        "galpop": galpop_file,
        "selfie": selfie_file,
        "catex": catex_file,
        "phase": phasedata_file,
        "tiles": tiles_file,
        "fibres": fibres_file
        }

    for _, items in file_dict.items():
        if items == '':
            raise ValueError('at least one sim file is missing')
        else: continue

    return file_dict

if __name__ == "__main__":
    import doctest
    doctest.testmod()