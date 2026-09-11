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


def adj_setup(data: Table, SNANA_temps:str) -> list:

    """
    Takes table of parameters extracted from simulations and
    adds then to lists for use in generating spectra.

    :param data: Table of parameters.
    :param SNANA_temps: path to SNANA SN SEDs
    :type data: Table
    :type SNANA_temps: str
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
    isky = []
    name = []
    ra = []
    dec = []
    x0s = []
    x1s = []
    cs = []
    model = []

    real_data = data

    for it in range(len(real_data)):
        Smags.append(real_data['mag'][it])
        Gmags.append(real_data['hostgal_mag_r'][it])
        redshift.append(real_data['redshift_estimate'][it])
        templates.append(real_data['template'][it])
        SNe_types.append(real_data['sim_type_name'][it])
        phase.append(-99)
        ddlr.append(real_data['hostgal_ddlr'][it])
        snsep.append(real_data['hostgal_snsep'][it])
        texp_visit.append(real_data['visit_texp'][it])
        isky.append(real_data['isky'][it])
        name.append(real_data['name'][it])
        ra.append(real_data['ra_y'][it])
        dec.append(real_data['dec_y'][it])
        x0s.append(real_data['sim_SALT2x0'][it])
        x1s.append(real_data['sim_SALT2x1'][it])
        cs.append(real_data['sim_SALT2c'][it])
        model.append(real_data['sim_type_name'][it])

    SNe_temp_array = glob.glob(SNANA_temps+'*')
    print(SNe_temp_array)

    SNe_temp_names = []
    for t in range(len(SNe_temp_array)):
        SNe_name_begins = SNe_temp_array[t].find('SNANA_temps/') + 12
        SNe_temp_names.append(SNe_temp_array[t][SNe_name_begins:])

    supernovae = []
    for ti in range(len(Smags)):

        supernovae.append(SNe_temp_array[SNe_temp_names.index(templates[ti])])



    print(len(Smags), len(Gmags), len(redshift), len(supernovae))

    return redshift, Smags, Gmags, phase, templates, SNe_types, supernovae, ddlr, snsep, texp_visit, isky, name, ra, dec, x0s, x1s, cs, model


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

def assign_host(trans_species, template_path):
    '''
    Assigns host galaxies to transients probabilistically based on the 
    probabilities within Hakobyan+2012. TDEs and CaRTs assigned random hosts
    and SLSNe always assigned Sc.

    :param trans_species: the subclass of the transient being operated on.
    :param template_path: path to the host template files.
    :type trans_species: str
    :type template_path: str
    :returns host: the path the host assigned for the input transient.
    "rtype: str
    '''

    host = 0

    # arrays here are the probability bins for a random seed to be a:
    # E, E/S0, S0, S0/a, Sa, Sab, Sb, Sbc, Sc, Scd, Sd+
    Ia = [0.08356164383561644, 0.1315068493150685, 0.20136986301369864, 0.3013698630136986, 0.34383561643835614, 0.41232876712328764, 0.5794520547945206, 0.7506849315068493, 0.9068493150684932, 0.9452054794520548, 1.0]
    Ib = [0.02040816326530612, 0.02040816326530612, 0.04081632653061224, 0.04081632653061224, 0.061224489795918366, 0.10204081632653061, 0.22448979591836735, 0.40816326530612246, 0.6938775510204082, 0.7959183673469388, 1.0]
    Ic = [0.0, 0.0, 0.0, 0.0, 0.011235955056179775, 0.056179775280898875, 0.2247191011235955, 0.5617977528089888, 0.7640449438202248, 0.8651685393258428, 1.0]
    II = [0.0, 0.0, 0.0037105751391465678, 0.011131725417439703, 0.027829313543599257, 0.04452690166975881, 0.22634508348794063, 0.45083487940630795, 0.7884972170686456, 0.8608534322820037, 1.0]
    IIn = [0.0, 0.0, 0.0, 0.013513513513513514, 0.04054054054054054, 0.04054054054054054, 0.1891891891891892, 0.3513513513513514, 0.7297297297297298, 0.8513513513513514, 1.0]
    IIb = [0.0, 0.0, 0.0, 0.0, 0.04878048780487805, 0.04878048780487805, 0.14634146341463417, 0.36585365853658536, 0.6097560975609756, 0.7317073170731707, 1.0]

    seed = random.random()

    if trans_species in ['Ia_', 'Iap']:
        #uses Ia probability distribution
        if seed < Ia[0]:
            # it's an ellipical
            host = glob.glob(str(template_path+'*el*'))
        elif Ia[0] <= seed < Ia[1]:
            #its either elliptical or S0
            host = random.choice([glob.glob(str(template_path+'*el*')),
                                  glob.glob(str(template_path+'*s0*'))])
        elif Ia[1] <= seed < Ia[2]:
            #its S0
            host = glob.glob(str(template_path+'*s0*'))
        elif Ia[2] <= seed < Ia[3]:
            #its either S0 or Sa
            host = random.choice([glob.glob(str(template_path+'*s0*')),
                                  glob.glob(str(template_path+'*sa*'))])    
        elif Ia[3] <= seed < Ia[4]:
            #its Sa
            host = glob.glob(str(template_path+'*sa*'))
        elif Ia[4] <= seed < Ia[5]:
            #its either Sa or Sb
            host = random.choice([glob.glob(str(template_path+'*sa*')),
                                  glob.glob(str(template_path+'*sb*'))])
        elif Ia[5] <= seed < Ia[6]:
            #its Sb
            host = glob.glob(str(template_path+'*sb*'))
        elif Ia[6] <= seed < Ia[7]:
            #its either Sb or Sc
            host = random.choice([glob.glob(str(template_path+'*sb*')),
                                  glob.glob(str(template_path+'*sc*'))])
            
        elif Ia[7] <= seed < Ia[8]:
            #its Sc
            host = glob.glob(str(template_path+'*sc*'))
        elif Ia[8] <= seed < Ia[9]:
            #its either Sc or Sd (we'll call this Sc)
            host = glob.glob(str(template_path+'*sc*'))

        elif seed > Ia[9]:
            #its Sd+ (so either Sc or Sb (the most common hosts))
            host = random.choice([glob.glob(str(template_path+'*sb*')),
                                  glob.glob(str(template_path+'*sc*'))])
        
        else:
            print(seed , ' is not recognised for some reason')

    elif trans_species == 'Ib_':
        #uses Ib probability distribution
        if seed < Ib[0]:
            # it's an ellipical
            host = glob.glob(str(template_path+'*el*'))
        elif Ib[0] <= seed < Ib[1]:
            #its either elliptical or S0
            host = random.choice([glob.glob(str(template_path+'*el*')),
                                  glob.glob(str(template_path+'*s0*'))])
        elif Ib[1] <= seed < Ib[2]:
            #its S0
            host = glob.glob(str(template_path+'*s0*'))
        elif Ib[2] <= seed < Ib[3]:
            #its either S0 or Sa
            host = random.choice([glob.glob(str(template_path+'*s0*')),
                                  glob.glob(str(template_path+'*sa*'))])    
        elif Ib[3] <= seed < Ib[4]:
            #its Sa
            host = glob.glob(str(template_path+'*sa*'))
        elif Ib[4] <= seed < Ib[5]:
            #its either Sa or Sb
            host = random.choice([glob.glob(str(template_path+'*sa*')),
                                  glob.glob(str(template_path+'*sb*'))])
        elif Ib[5] <= seed < Ib[6]:
            #its Sb
            host = glob.glob(str(template_path+'*sb*'))
        elif Ib[6] <= seed < Ib[7]:
            #its either Sb or Sc
            host = random.choice([glob.glob(str(template_path+'*sb*')),
                                  glob.glob(str(template_path+'*sc*'))])
            
        elif Ib[7] <= seed < Ib[8]:
            #its Sc
            host = glob.glob(str(template_path+'*sc*'))
        elif Ib[8] <= seed < Ib[9]:
            #its either Sc or Sd (we'll call this Sc)
            host = glob.glob(str(template_path+'*sc*'))

        elif seed > Ib[9]:
            #its Sd+ (so either Sc or Sb (the most common hosts))
            host = random.choice([glob.glob(str(template_path+'*sb*')),
                                  glob.glob(str(template_path+'*sc*'))])
        
        else:
            print(seed , ' is not recognised for some reason')

    elif trans_species == 'Ic_':
        #uses Ic probability distribution
        if seed < Ic[0]:
            # it's an ellipical
            host = glob.glob(str(template_path+'*el*'))
        elif Ic[0] <= seed < Ic[1]:
            #its either elliptical or S0
            host = random.choice([glob.glob(str(template_path+'*el*')),
                                  glob.glob(str(template_path+'*s0*'))])
        elif Ic[1] <= seed < Ic[2]:
            #its S0
            host = glob.glob(str(template_path+'*s0*'))
        elif Ic[2] <= seed < Ic[3]:
            #its either S0 or Sa
            host = random.choice([glob.glob(str(template_path+'*s0*')),
                                  glob.glob(str(template_path+'*sa*'))])    
        elif Ic[3] <= seed < Ic[4]:
            #its Sa
            host = glob.glob(str(template_path+'*sa*'))
        elif Ic[4] <= seed < Ic[5]:
            #its either Sa or Sb
            host = random.choice([glob.glob(str(template_path+'*sa*')),
                                  glob.glob(str(template_path+'*sb*'))])
        elif Ic[5] <= seed < Ic[6]:
            #its Sb
            host = glob.glob(str(template_path+'*sb*'))
        elif Ic[6] <= seed < Ic[7]:
            #its either Sb or Sc
            host = random.choice([glob.glob(str(template_path+'*sb*')),
                                  glob.glob(str(template_path+'*sc*'))])
            
        elif Ic[7] <= seed < Ic[8]:
            #its Sc
            host = glob.glob(str(template_path+'*sc*'))
        elif Ic[8] <= seed < Ic[9]:
            #its either Sc or Sd (we'll call this Sc)
            host = glob.glob(str(template_path+'*sc*'))

        elif seed > Ic[9]:
            #its Sd+ (so either Sc or Sb (the most common hosts))
            host = random.choice([glob.glob(str(template_path+'*sc*')),
                                  glob.glob(str(template_path+'*sb*'))])
        
        else:
            print(seed , ' is not recognised for some reason')

    elif trans_species == 'II_':
        print(II[2], II[3])
        #uses II probability distribution
        if seed < II[0]:
            # it's an ellipical
            host = glob.glob(str(template_path+'*el*'))
        elif II[0] <= seed < II[1]:
            #its either elliptical or S0
            host = random.choice([glob.glob(str(template_path+'*el*')),
                                  glob.glob(str(template_path+'*s0*'))])
        elif II[1] <= seed < II[2]:
            #its S0
            host = glob.glob(str(template_path+'*s0*'))
        elif II[2] <= seed < II[3]:
            #its either S0 or Sa
            host = random.choice([glob.glob(str(template_path+'*s0*')),
                                  glob.glob(str(template_path+'*sa*'))])    
        elif II[3] <= seed < II[4]:
            #its Sa
            host = glob.glob(str(template_path+'*sa*'))
        elif II[4] <= seed < II[5]:
            #its either Sa or Sb
            host = random.choice([glob.glob(str(template_path+'*sa*')),
                                  glob.glob(str(template_path+'*sb*'))])
        elif II[5] <= seed < II[6]:
            #its Sb
            host = glob.glob(str(template_path+'*sb*'))
        elif II[6] <= seed < II[7]:
            #its either Sb or Sc
            host = random.choice([glob.glob(str(template_path+'*sb*')),
                                  glob.glob(str(template_path+'*sc*'))])
            
        elif II[7] <= seed < II[8]:
            #its Sc
            host = glob.glob(str(template_path+'*sc*'))
        elif II[8] <= seed < II[9]:
            #its either Sc or Sd (we'll call this Sc)
            host = glob.glob(str(template_path+'*sc*'))

        elif seed > II[9]:
            #its Sd+ (so either Sc or Sb (the most common hosts))
            host = random.choice([glob.glob(str(template_path+'*sb*')),
                                  glob.glob(str(template_path+'*sc*'))])
        
        else:
            print(seed , ' is not recognised for some reason')

    elif trans_species == 'IIn':
        #uses IIn probability distribution
        if seed < IIn[0]:
            # it's an ellipical
            host = glob.glob(str(template_path+'*el*'))
        elif IIn[0] <= seed < IIn[1]:
            #its either elliptical or S0
            host = random.choice([glob.glob(str(template_path+'*el*')),
                                  glob.glob(str(template_path+'*s0*'))])
        elif IIn[1] <= seed < IIn[2]:
            #its S0
            host = glob.glob(str(template_path+'*s0*'))
        elif IIn[2] <= seed < IIn[3]:
            #its either S0 or Sa
            host = random.choice([glob.glob(str(template_path+'*s0*')),
                                  glob.glob(str(template_path+'*sa*'))])    
        elif IIn[3] <= seed < IIn[4]:
            #its Sa
            host = glob.glob(str(template_path+'*sa*'))
        elif IIn[4] <= seed < IIn[5]:
            #its either Sa or Sb
            host = random.choice([glob.glob(str(template_path+'*sa*')),
                                  glob.glob(str(template_path+'*sb*'))])
        elif IIn[5] <= seed < IIn[6]:
            #its Sb
            host = glob.glob(str(template_path+'*sb*'))
        elif IIn[6] <= seed < IIn[7]:
            #its either Sb or Sc
            host = random.choice([glob.glob(str(template_path+'*sb*')),
                                  glob.glob(str(template_path+'*sc*'))])
            
        elif IIn[7] <= seed < IIn[8]:
            #its Sc
            host = glob.glob(str(template_path+'*sc*'))
        elif IIn[8] <= seed < IIn[9]:
            #its either Sc or Sd (we'll call this Sc)
            host = glob.glob(str(template_path+'*sc*'))

        elif seed > IIn[9]:
            #its Sd+ (so either Sc or Sb (the most common hosts))
            host = random.choice([glob.glob(str(template_path+'*sc*')),
                                  glob.glob(str(template_path+'*sb*'))])
        
        else:
            print(seed , ' is not recognised for some reason')

    elif trans_species == 'IIb':
        #uses IIb probability distribution
        if seed < IIb[0]:
            # it's an ellipical
            host = glob.glob(str(template_path+'*el*'))
        elif IIb[0] <= seed < IIb[1]:
            #its either elliptical or S0
            host = random.choice([glob.glob(str(template_path+'*el*')),
                                  glob.glob(str(template_path+'*s0*'))])
        elif IIb[1] <= seed < IIb[2]:
            #its S0
            host = glob.glob(str(template_path+'*s0*'))
        elif IIb[2] <= seed < IIb[3]:
            #its either S0 or Sa
            host = random.choice([glob.glob(str(template_path+'*s0*')),
                                  glob.glob(str(template_path+'*sa*'))])    
        elif IIb[3] <= seed < IIb[4]:
            #its Sa
            host = glob.glob(str(template_path+'*sa*'))
        elif IIb[4] <= seed < IIb[5]:
            #its either Sa or Sb
            host = random.choice([glob.glob(str(template_path+'*sa*')),
                                  glob.glob(str(template_path+'*sb*'))])
        elif IIb[5] <= seed < IIb[6]:
            #its Sb
            host = glob.glob(str(template_path+'*sb*'))
        elif IIb[6] <= seed < IIb[7]:
            #its either Sb or Sc
            host = random.choice([glob.glob(str(template_path+'*sb*')),
                                  glob.glob(str(template_path+'*sc*'))])
            
        elif IIb[7] <= seed < IIb[8]:
            #its Sc
            host = glob.glob(str(template_path+'*sc*'))
        elif IIb[8] <= seed < IIb[9]:
            #its either Sc or Sd (we'll call this Sc)
            host = glob.glob(str(template_path+'*sc*'))

        elif seed > IIb[9]:
            #its Sd+ (so either Sc or Sb (the most common hosts))
            host = random.choice([glob.glob(str(template_path+'*sc*')),
                                  glob.glob(str(template_path+'*sb*'))])
        
        else:
            print(seed , ' is not recognised for some reason')

    elif trans_species == 'SL_':
        #always an Sc
        host = glob.glob(str(template_path+'*sc*'))

    elif trans_species in ['TDE', 'CRT', 'KN_']:
        #just random
        host = [random.choice(glob.glob(str(template_path+'*')))]

    else:
        print(trans_species, ' species isnt recognised for some reason')

    return host

def get_SN_and_host_type(gal_mags, SN, model_name, host_location):

    SN_type_str = []
    SN_phase = []
    galaxies = []

    for j in range(len(gal_mags)):
        SN_ID_str_begin = SN[j].find('snt') + 3
        SN_ID_str_end = SN[j].find('_phase')
        SN_phase_end = SN[j].find('_redshift')
        snt_ID = int(SN[j][SN_ID_str_begin:SN_ID_str_end])
        SN_phase.append(int(SN[j][SN_ID_str_end+6:SN_phase_end]))

        if snt_ID == 1:
        #for some reason at least one IIb has a spec with SNT==1 -> Ia???? breaks SALT2 generation
            if model_name[j] == 'SALT2.WF':
                SN_type_str.append('Ia_')
                SN[j] = 'SALT2'
            else:
                SN_type_str.append('IIb')
        elif snt_ID in [11, 12]:
            SN_type_str.append('Iap')
        elif snt_ID == 60:
            SN_type_str.append('KN_')
        elif snt_ID == 70:
            SN_type_str.append('SL_')
        elif snt_ID == 50:
            SN_type_str.append('CRT')
        elif snt_ID == 80:
            SN_type_str.append('TDE')
        elif snt_ID == 21:
            SN_type_str.append('IIn')
        elif snt_ID == 23:
            SN_type_str.append('IIb')
        elif snt_ID == 25:
            SN_type_str.append('II_')
        elif snt_ID == 32:
            SN_type_str.append('Ib_')
        elif snt_ID in [33, 35]:
            SN_type_str.append('Ic_')
        elif snt_ID == 20:
            SN_type_str.append('CC_')
        else:
            print(f'Unknown snt ID string: {snt_ID}, I would look into this')

        galaxies.append(assign_host(SN_type_str[-1], host_location)[0])

        return galaxies, SN_type_str, SN_phase

if __name__ == "__main__":
    import doctest
    doctest.testmod()