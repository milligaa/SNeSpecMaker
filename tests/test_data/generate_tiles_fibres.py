from astropy.io import fits
from astropy.table import Table
import numpy as np

hdul = fits.open('/Users/andrew/Desktop/Python_Stuff/my_spectrum_maker/SNeSpecMaker/tests/test_data/S238_fibres_test.fits')
print(hdul[1].columns)

# targ_ids = [50371430,50374278,50376190]

# good_index = []
# targ_id_list = list(hdul[1].data['targ_id'])

# for i in range(len(targ_ids)):
#     good_index.append(targ_id_list.index(targ_ids[i]))

# print(hdul[1].data['tile_id'][good_index[0]])
# print(hdul[1].data['tile_id'][good_index[1]])
# print(hdul[1].data['tile_id'][good_index[2]])

# to_save = Table()
# to_save['tile_id'] = [hdul[1].data['tile_id'][good_index[0]],
#                       hdul[1].data['tile_id'][good_index[1]],
#                       hdul[1].data['tile_id'][good_index[2]]]
# to_save['targ_id'] = targ_ids

# to_save.write('/Users/andrew/Desktop/Python_Stuff/my_spectrum_maker/SNeSpecMaker/tests/test_data/S238_fibres_test.fits', format='fits')