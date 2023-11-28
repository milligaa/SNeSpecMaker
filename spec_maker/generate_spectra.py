from astropy.table import Table
from SNeSpecMaker.Utils.sim_observation import Comb_Maker
from SNeSpecMaker.Utils.param_setup import adj_setup, SELFIE_extractor
from SNeSpecMaker.Utils.seeing_effects import point_convolute
from SNeSpecMaker.Utils.seeing_effects import effective_fibre_mag
import argparse
import yaml


def template_flow(yaml_init_data: dict) -> None:
    """
    This functions combines those in the util folder to create
    contaminated spectra and save the SNR data from them.

    :param str yaml_init_data: dictionary of paths loaded from yaml file
    """

    results_path = yaml_init_data["SNR_results_path"]
    sim_data_location = yaml_init_data["sim_data_loc"]
    spectrum_save_path = yaml_init_data["spectra_saving_path"]
    save_tab_path = yaml_init_data["interim_product_path"]

    # extracts data from the selfie files and returns
    # the table with everything I need

    data_to_use = SELFIE_extractor(sim_data_location, save_tab_path)

    # uses the SELFIE file to get parameters like redshift,
    # magnitude and SNe type for use in template generation
    template_values = adj_setup(data_to_use)
    redshift = template_values[0]
    Smags = template_values[1]
    Gmags = template_values[2]
    phase = template_values[3]
    galaxies = template_values[6]
    supernovae = template_values[7]
    ddlr = template_values[8]
    snsep = template_values[9]

    # now must perform the correction for effective fibre mag
    seeing_val = 0.8
    gmag_eff_fibre = []
    smag_eff_fibre = []
    for h in range(len(Gmags)):
        gmag_eff_fibre.append(
            effective_fibre_mag(
                snsep[h], ddlr[h], 0.5, Gmags[h], snsep[h] / 100, seeing_val
            )
        )

        smag_eff_fibre.append(point_convolute(seeing=seeing_val,
                                              sne_mag=Smags[h]))

    for o in range(len(gmag_eff_fibre)):
        Gmags[o] = gmag_eff_fibre[o]
        Smags[o] = smag_eff_fibre[o]

    L1_SNR_corr1 = []
    Comb_mag_corr1 = []

    L1_SNR_corr2 = []
    Comb_mag_corr2 = []

    L1_SNR_corr3 = []
    Comb_mag_corr3 = []

    L1_SNR_corr4 = []
    Comb_mag_corr4 = []

    bad_index = []
    SNR_append_index = []

    SN_type_str = []

    for j in range(len(Gmags)):
        SN_ID_str_begin = supernovae[j].find("snt") + 3
        SN_ID_str_end = supernovae[j].find("_phase")
        snt_ID = int(supernovae[j][SN_ID_str_begin:SN_ID_str_end])

        if snt_ID in [1, 11, 12]:
            SN_type_str.append("Ia_")
        elif snt_ID == 60:
            SN_type_str.append("KN_")
        elif snt_ID == 70:
            SN_type_str.append("SL_")
        elif snt_ID == 50:
            SN_type_str.append("CRT")
        elif snt_ID == 80:
            SN_type_str.append("TDE")
        elif snt_ID == 21:
            SN_type_str.append("IIn")
        elif snt_ID == 23:
            SN_type_str.append("IIb")
        elif snt_ID == 25:
            SN_type_str.append("II_")
        elif snt_ID == 32:
            SN_type_str.append("Ib_")
        elif snt_ID in [33, 35]:
            SN_type_str.append("Ic_")
        elif snt_ID == 20:
            SN_type_str.append("CC_")
        else:
            print("damn")

    for dummy in range(len(Gmags)):
        if len(L1_SNR_corr1) == len(L1_SNR_corr4):
            try:
                print(dummy, "of", len(supernovae))
                result = Comb_Maker(
                    supernovae[dummy],
                    galaxies[dummy],
                    50,
                    Smags[dummy],
                    redshift[dummy],
                    galaxies[dummy][72:74],
                    SN_type_str[dummy],
                    20,
                    seeing_val,
                    spectrum_save_path,
                )

                result2 = Comb_Maker(
                    supernovae[dummy],
                    galaxies[dummy],
                    Gmags[dummy],
                    Smags[dummy],
                    redshift[dummy],
                    galaxies[dummy][72:74],
                    SN_type_str[dummy],
                    40,
                    seeing_val,
                    spectrum_save_path,
                )

                result3 = Comb_Maker(
                    supernovae[dummy],
                    galaxies[dummy],
                    50,
                    Smags[dummy],
                    redshift[dummy],
                    galaxies[dummy][72:74],
                    SN_type_str[dummy],
                    60,
                    seeing_val,
                    spectrum_save_path,
                )

                result4 = Comb_Maker(
                    supernovae[dummy],
                    galaxies[dummy],
                    Gmags[dummy],
                    Smags[dummy],
                    redshift[dummy],
                    galaxies[dummy][72:74],
                    SN_type_str[dummy],
                    120,
                    seeing_val,
                    spectrum_save_path,
                )

                L1_SNR_corr1.append(result[0])
                Comb_mag_corr1.append(result[1].value)

                L1_SNR_corr2.append(result2[0])
                Comb_mag_corr2.append(result2[1].value)

                L1_SNR_corr3.append(result3[0])
                Comb_mag_corr3.append(result3[1].value)

                L1_SNR_corr4.append(result4[0])
                Comb_mag_corr4.append(result4[1].value)

                SNR_append_index.append(dummy)

            except Exception as excpt:
                print("duplicate", dummy)
                bad_index.append(dummy)
                print(excpt)

        else:
            print("array lenghts diverge here", dummy)
            print(
                len(L1_SNR_corr1),
                len(L1_SNR_corr2),
                len(L1_SNR_corr3),
                len(L1_SNR_corr4),
            )
            print(result, result2, result3, result4)
            break

    for bad in range(len(bad_index)):
        del redshift[bad_index[bad] - bad]
        del Smags[bad_index[bad] - bad]
        del Gmags[bad_index[bad] - bad]
        del phase[bad_index[bad] - bad]

    SNR_table = Table()
    SNR_table["average_SNR"] = L1_SNR_corr1
    SNR_table["combined_mag"] = Comb_mag_corr1
    SNR_table["Redshift"] = redshift
    SNR_table["Smags"] = Smags
    SNR_table["gmags"] = Gmags
    SNR_table["phase"] = phase

    SNR_table2 = Table()
    SNR_table2["average_SNR"] = L1_SNR_corr2
    SNR_table2["combined_mag"] = Comb_mag_corr2
    SNR_table2["Redshift"] = redshift
    SNR_table2["Smags"] = Smags
    SNR_table2["gmags"] = Gmags
    SNR_table2["phase"] = phase

    SNR_table3 = Table()
    SNR_table3["average_SNR"] = L1_SNR_corr3
    SNR_table3["combined_mag"] = Comb_mag_corr3
    SNR_table3["Redshift"] = redshift
    SNR_table3["Smags"] = Smags
    SNR_table3["gmags"] = Gmags
    SNR_table3["phase"] = phase

    SNR_table4 = Table()
    SNR_table4["average_SNR"] = L1_SNR_corr4
    SNR_table4["combined_mag"] = Comb_mag_corr4
    SNR_table4["Redshift"] = redshift
    SNR_table4["Smags"] = Smags
    SNR_table4["gmags"] = Gmags
    SNR_table4["phase"] = phase

    ascii.write(
        SNR_table,
        (results_path + "S238missed_wfd_texpobj_nohost.txt"),
        overwrite=True
    )
    ascii.write(
        SNR_table2,
        (results_path + "S238missed_wfd_texpobj_host.txt"),
        overwrite=True
    )
    ascii.write(
        SNR_table3,
        (results_path + "S238missed_wfd_texpvis_nohost.txt"),
        overwrite=True
    )
    ascii.write(
        SNR_table4,
        (results_path + "S238missed_wfd_texpvis_host.txt"),
        overwrite=True
    )


def parser():
    """
    Function that allows a filename to be passed from command line
    into the template_flow function

    :return: object containing arguments passed in the command line as strings
    :rtype: argparse.ArgmumentParser class
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--initfile",
        default=None,
        help="yaml file containing input frames and other procedure arguments",
    )
    return parser


if __name__ == "__main__":
    parser = parser()
    args = parser.parse_args()
    print(args.initfile)
    with open(str(args.initfile), "r") as file:
        params = yaml.safe_load(file)
