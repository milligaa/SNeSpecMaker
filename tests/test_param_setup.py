import pytest
from ..Utils.param_setup import assign_sim_files

good_data_path = 'good_test_data/'

# tests of assign_sim_files
def test_with_good_data():
    output = assign_sim_files(good_data_path)
    assert output["galpop"] == str(good_data_path+'3pt2_wfd_host_test.csv')
    assert output["tiles"] == str(good_data_path+'S238_tiles_test.fits')
                                   