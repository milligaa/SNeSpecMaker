#when called will remove the json file called parameters.json so that the next parameter file can then be renamed to parameters.json

import os
import glob
import json

#gtet all jsons left
jsons = glob.glob('/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/NGSF-multi/config/*.json')
jsons.sort()

#remove the one just used
os.remove(jsons[0])


#get data and make a new parameter file called parameters.json and then delete the source
with open(jsons[1], 'r') as file:
    data = json.load(file)

with open('/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/NGSF-multi/config/parameters.json', 'w') as file:
    json.dump(data, file, indent=4)


os.remove(jsons[1])
