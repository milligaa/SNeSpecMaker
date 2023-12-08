import sys
from NGSF.sf_class import Superfit
import json

#load json file and get object file str from it
with open('/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/NGSF-multi/config/parameters.json', 'r') as file:
    data = json.load(file)
    file_str = str(data['object_to_fit'])

supernova = Superfit(file_str)
supernova.superfit()
