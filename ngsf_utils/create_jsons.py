#this code finds all spectra wanting to be fit and then creates a json file with correct filename and redshift range (small range around true redshift)

import glob
import json

print('running create_many.py')

spectra = glob.glob('/Users/andrew/Desktop/Python_Stuff/Contaminated_Templates/S238missed_expvisit_host/*')
spectra.sort()

print('there are, ', len(spectra), ' spectra')

#create list of redshifts
redshifts = []
for i in range(len(spectra)):
    start = spectra[i].find('z') + 1
    end = spectra[i].find('texp')
    redshifts.append(spectra[i][start:end])

print(redshifts, ' is the redshifts')


#load parameters file data
for j in range(len(spectra)):
    with open('/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/NGSF-multi/config/parameters.json', 'r') as file:
        data = json.load(file)
        data["object_to_fit"] = "{0}".format(spectra[j])
        data["use_exact_z"] = 0
        data["z_exact"] = float(redshifts[j])
        data["z_range_begin"] = round(float(redshifts[j]) -0.01, 6)
        data["z_range_end"] = round(float(redshifts[j]) +0.01, 6)
        data['z_int'] = 0.01
        data["error_spectrum"] = "attached"
        data["lower_lam"] = 4000
        data["upper_lam"] = 8500

    new_json_str = str('/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/NGSF-multi/config/parameters' + str(j) + '.json')

    with open(new_json_str, 'w') as f:
        json.dump(data, f, indent=4)

print('run_many.py complete')