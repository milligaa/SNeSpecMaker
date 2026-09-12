import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

from SNeSpecMaker.Utils.sim_observation import Comb_Maker
from SNeSpecMaker.Utils.param_setup import adj_setup, get_SN_and_host_type
from SNeSpecMaker.Utils.seeing_effects import point_convolute
from SNeSpecMaker.Utils.seeing_effects import effective_fibre_mag
from astropy.table import Table, join
import argparse
import yaml
from astropy.io import ascii


def make_blended(save_spec_loc, population_file, begin, end, host_loc,
                 SNANA_temp_path):

    data_to_use = Table.read(population_file,
                             format = 'csv', delimiter = ',')[begin:end]

    print('----------------------------------')
    print(f'Number of spectra: {len(data_to_use)}')
    print('----------------------------------')

    #uses the SELFIE file to get parameters like redshift, magnitude and SNe type for use in template generation
    template_values = adj_setup(data_to_use, SNANA_temp_path)
    redshift = template_values[0]
    Smags = template_values[1]
    Gmags = template_values[2]
    supernovae = template_values[5]
    ddlr = template_values[6]
    snsep = template_values[7]
    texp_visit = template_values[8]
    name = template_values[9]
    ra = template_values[10]
    dec = template_values[11]
    x0 = template_values[12]
    x1 = template_values[13]
    c = template_values[14]
    models = template_values[15]

    #now must perform the correction for effective fibre mag
    seeing_val = 0.8
    gmag_eff_fibre = []
    smag_eff_fibre = []
    for h in range(len(Gmags)):
        print('performing magnitude correction ', h+1, ' of ', len(Gmags))
        gmag_eff_fibre.append(effective_fibre_mag(snsep[h], ddlr[h], 0.5, Gmags[h], snsep[h]/100, seeing_val))
        smag_eff_fibre.append(point_convolute(seeing=seeing_val, sne_mag=Smags[h]))

    L1_SNR_corr1 = []
    Comb_mag_corr1 = []

    SNR_append_index = []

    galaxies, SN_type_str, SN_phase = get_SN_and_host_type(Gmags, supernovae,
                                                           models, host_loc)

    gal_type_str = []
    for dummy in range(len(Gmags)):

        print(supernovae[dummy])

        #construct salt2 param dictionary
        salt2_pars = {'x0': x0[dummy], 'x1':x1[dummy],
                      'c': c[dummy], 'phase': SN_phase[dummy]}
            
        gal_str_start = galaxies[dummy].find('kinney/')
        gal_str_end = galaxies[dummy].find('_template.fits')
        gal_type = galaxies[dummy][gal_str_start+7:gal_str_end]
        gal_type_str.append(gal_type)
        print(dummy, 'of', len(supernovae))
        result = Comb_Maker(supernovae[dummy], galaxies[dummy],
        gmag_eff_fibre[dummy], smag_eff_fibre[dummy], redshift[dummy], gal_type, SN_type_str[dummy],
        texp_visit[dummy], name[dummy], seeing_val,
        save_spec_loc, SALT2_params=salt2_pars)

        L1_SNR_corr1.append(result[0])
        Comb_mag_corr1.append(result[1].value)

        SNR_append_index.append(dummy)
                                                                               

    SNR_table = Table()
    SNR_table['Combined_SNR'] = L1_SNR_corr1
    SNR_table['Combined_mag'] = Comb_mag_corr1
    SNR_table['Redshift'] = redshift
    SNR_table['True_Smag'] = Smags
    SNR_table['True_Gmag'] = Gmags
    SNR_table['Fibre_Smag'] = smag_eff_fibre
    SNR_table['Fibre_Gmag'] = gmag_eff_fibre
    SNR_table['Phase'] = SN_phase
    SNR_table['Texp_Visit'] = texp_visit
    SNR_table['sne_class'] = gal_type_str
    SNR_table['host_morph'] = SN_type_str
    SNR_table['name'] = name
    SNR_table['ra'] = ra
    SNR_table['dec'] = dec

    ascii.write(SNR_table, str(save_spec_loc+'blended_res'+str(begin)+'.csv'), overwrite = True)

def parser():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--initfile",
        default=None,
        help="location of initialisation file (yaml)",
    )

    return parser

if __name__ == "__main__":
    parser = parser()
    args = parser.parse_args()
    print(args.initfile)
    with open(str(args.initfile), "r") as file:
        params = yaml.safe_load(file)

    make_blended(params['spectra_save_path'], params['input_population'],
                 params['begin'], params['end'], params['host_loc'],
                 params['SNANA_SED_loc'])