import time
import glob
import astrodash
import numpy as np
from astropy.table import Table
from astropy.io import ascii

def big_DASH(filename, redshift, results_path):
    print('starting big DASH for ', filename, ' exposures')

    #first step is to populate an array with the filenames and redshifts


    start = time.time()

    glob_search_name = str(filename) + '/*'

    files = glob.glob(glob_search_name)
    files.sort()

    #this gets the redshifts from the filenames
    redshifts = []
    Smag = []
    Gmag = []
    SNe = []
    Gtype = []
    texp = []

    for i in range(len(files)):
        gmag = files[i].find('Gmag')
        smag = files[i].find('Smag')
        z = files[i].find('z')
        texp_str = files[i].find('texp')
        txt_str = files[i].find('.txt')
        redshifts.append(float(files[i][z+1:texp_str]))
        Smag.append(float(files[i][smag+4:gmag]))
        Gmag.append(float(files[i][gmag+4:z]))
        SNe.append(files[i][smag-3:smag])
        Gtype.append(files[i][smag-5:smag-3])
        texp.append(files[i][texp_str+4:txt_str])

    example = []

#creates array of files and corresponding redhift in format required for dash
    for x in range(len(files)):
        example.append((files[x], redshifts[x]))

    # Create filenames and knownRedshifts lists
    filenames = [i[0] for i in example]
    knownRedshifts = [i[1] for i in example]

    # Classify all spectra
    classification = astrodash.Classify(filenames, knownRedshifts, classifyHost=False, knownZ=redshift, smooth=6, rlapScores = True)
    bestFits, Redshifts, bestTypes, rlapFlag, matchesFlag, redshiftErrs = classification.list_best_matches(n=5)

    comb_table = Table()
    comb_table['SNe'] = SNe
    comb_table['Smag'] = Smag
    comb_table['Gmag'] = Gmag
    comb_table['BestFit'] = bestFits
    comb_table['BestType'] = bestTypes
    comb_table['mathflag'] = matchesFlag
    comb_table['Gal'] = Gtype
    comb_table['z'] = redshifts
    comb_table['rlaps'] = rlapFlag
    comb_table['texp'] = texp

    if redshift == True:
        filename_forsave = str(results_path)+'/DASH_'+ str(filename) +'knownz.txt'
    elif redshift == False:
        filename_forsave = str(results_path)+'/DASH_'+ str(filename) +'unknownz.txt'
    # filename_forsave = '/Users/Andrew/Desktop/Python_Stuff/SN_and_Galaxy/results/text_file_dump4/big_DASH_SEFLIE172_allres_noz.txt'
    ascii.write(comb_table, filename_forsave, overwrite = True, format = 'ecsv')
    # Plot sn2013fs from open supernova catalog (2nd spectrum)
    # classification.plot_with_gui(indexToPlot=0)

    end = time.time()
    print('time for ', len(files), ' objects is', end - start, ' seconds')