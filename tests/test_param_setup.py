import pytest
from ..Utils.param_setup import assign_sim_files

good_data_path = 'param_setup_test_data/good_test_data/'
missing_file_path = 'param_setup_test_data/test_data_missing_file/'
wrong_filetype_path = 'param_setup_test_data/wrong_filetype.txt'
csv_wrong_header_path = 'param_setup_test_data/csv_wrong_header.csv'
fits_wrong_header_path = 'param_setup_test_data/fits_wrong_header.fits'
non_text_data_path = 'param_setup_test_data/test_non_text_data.png'

# tests of assign_sim_files
def test_with_good_data():
    output = assign_sim_files(good_data_path)
    assert output["galpop"] == str(good_data_path+'3pt2_wfd_host_test.csv')
    assert output["tiles"] == str(good_data_path+'S238_tiles_test.fits')

def test_with_missing_sim_file():
    with pytest.raises(ValueError):
        assign_sim_files(missing_file_path)

def test_with_text_file():
    with pytest.raises(IndexError):
        assign_sim_files(wrong_filetype_path)

def test_csv_wrong_header():
    with pytest.raises(IndexError):
        assign_sim_files(csv_wrong_header_path)

def test_wrong_fits_file():
    with pytest.raises(IndexError):
        assign_sim_files(fits_wrong_header_path)

def test_with_non_textbased_data():
    with pytest.raises(TypeError):
        assign_sim_files(non_text_data_path)

# testing the selfie extractor
