#this program takes the results from superfit and concatenates the best result from each into a single csv file

from astropy.io import ascii
from astropy.table import Table
import glob
import argparse

def ngsf_concat(csv_loc, result_path):
    #load all the csv files
    csv_files = glob.glob(str(csv_loc+'*.csv'))
    print(len(csv_files))

    final_res_table = Table()
    spectrum = []
    galfit_1 = []
    SNfit_1 = []
    z_1 = []
    redchi_1 = []
    SN_frac_1 = []
    phase_1 = []
    Av_1 = []
    galfit_2 = []
    SNfit_2 = []
    z_2 = []
    redchi_2 = []
    SN_frac_2 = []
    phase_2 = []
    Av_2 = []
    galfit_3 = []
    SNfit_3 = []
    z_3 = []
    redchi_3 = []
    SN_frac_3 = []
    phase_3 = []
    Av_3 = []
    texp = []
    redredchi_1 = []

    for i in range(len(csv_files)):
        print(i)
        table = Table.read(csv_files[i], format = 'csv', delimiter = ',')
        
        best_fit = table[0]

        spectrum.append(best_fit[0])

        start = best_fit[0].find('texp') + 4
        end = best_fit[0].find('host')

        texp.append(best_fit[0][start:end])
        galfit_1.append(best_fit[1])
        z_1.append(best_fit[5])
        phase_1.append(best_fit[7])
        Av_1.append(best_fit[6])
        redchi_1.append(best_fit[11])
        redredchi_1.append(best_fit[12])
        SN_frac_1.append(best_fit[9])

        SNefit_str_end = best_fit[2].find('phase')
        SNfit_1.append(best_fit[2][:SNefit_str_end])

        try:
            second_fit = table[1]

            galfit_2.append(second_fit[1])
            z_2.append(second_fit[5])
            phase_2.append(second_fit[7])
            Av_2.append(second_fit[6])
            redchi_2.append(second_fit[11])
            SN_frac_2.append(second_fit[9])

            SNefit_str_end = second_fit[2].find('phase')
            SNfit_2.append(second_fit[2][:SNefit_str_end])
        except:
            galfit_2.append('none')
            z_2.append('none')
            phase_2.append('none')
            Av_2.append('none')
            redchi_2.append('none')
            SN_frac_2.append('none')
            SNfit_2.append('none')

        try:
            third_fit = table[2]

            galfit_3.append(third_fit[1])
            z_3.append(third_fit[5])
            phase_3.append(third_fit[7])
            Av_3.append(third_fit[6])
            redchi_3.append(third_fit[11])
            SN_frac_3.append(third_fit[9])

            SNefit_str_end = third_fit[2].find('phase')
            SNfit_3.append(third_fit[2][:SNefit_str_end])    
        except:
            galfit_3.append('none')
            z_3.append('none')
            phase_3.append('none')
            Av_3.append('none')
            redchi_3.append('none')
            SN_frac_3.append('none')
            SNfit_3.append('none')

        print(spectrum[i], z_1[i])

    variables = [spectrum,
    galfit_1,
    SNfit_1,
    z_1,
    redchi_1,
    redredchi_1,
    SN_frac_1,
    phase_1,
    Av_1,
    galfit_2,
    SNfit_2,
    z_2,
    redchi_2,
    SN_frac_2,
    phase_2,
    Av_2,
    galfit_3,
    SNfit_3,
    z_3,
    redchi_3,
    SN_frac_3,
    phase_3,
    Av_3,
    texp]

    var_names = ['spectrum',
    'galfit_1',
    'SNfit_1',
    'z_1',
    'redchi_1',
    'redredchi_1',
    'SN_frac_1',
    'phase_1',
    'Av_1',
    'galfit_2',
    'SNfit_2',
    'z_2',
    'redchi_2',
    'SN_frac_2',
    'phase_2',
    'Av_2',
    'galfit_3',
    'SNfit_3',
    'z_3',
    'redchi_3',
    'SN_frac_3',
    'phase_3',
    'Av_3',
    'texp']

    for j in range(len(variables)):
        final_res_table[str(var_names[j])] = variables[j]

    ascii.write(final_res_table, result_path, 
                format = 'csv', delimiter = ',', overwrite=True)
    
def parser():
    """
    Function that allows a filename to be passed from command line
    into the template_flow function

    :return: object containing arguments passed in the command line as strings
    :rtype: argparse.ArgmumentParser class
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-c",
        "--csv",
        default=None,
        help="location of ngsf classification results",
    )

    parser.add_argument(
        "-r",
        "--results",
        default=None,
        help="path to save concatenated results",
    )

    return parser


if __name__ == "__main__":
    parser = parser()
    args = parser.parse_args()
    print(args.initfile)

