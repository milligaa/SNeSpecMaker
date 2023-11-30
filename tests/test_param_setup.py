import pytest
from ..Utils.param_setup import assign_sim_files, SELFIE_extractor, adj_setup

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
@pytest.fixture
def assign_good_sim_files():
    return assign_sim_files(good_data_path)


def test_selfie_extract_good_data(assign_good_sim_files):
    table = SELFIE_extractor(
        assign_good_sim_files,
        '/Users/andrew/Desktop/Python_Stuff/my_spectrum_maker/test_folder/combined_tables'
        )
    
    assert len(table.columns) == 31
    assert len(table) == 3
    assert table[1]['redshift_estimate'] == 0.545628


def test_extract_empty_dictionary_csv(assign_good_sim_files):
    
    assign_good_sim_files["galpop"] = ''

    with pytest.raises(FileNotFoundError):
        SELFIE_extractor(
            assign_good_sim_files,
            'test_folder/SNR_res/'
            )


def test_extract_empty_dictionary_fits(assign_good_sim_files):
    
    assign_good_sim_files["fibres"] = ''

    with pytest.raises(FileNotFoundError):
        SELFIE_extractor(
            assign_good_sim_files,
            'test_folder/SNR_res/'
            )

def test_extract_csv_wrong_columns(assign_good_sim_files):

    assign_good_sim_files["selfie"] = 'param_setup_test_data/S238_S1001_test_missing_col.csv'

    with pytest.raises(KeyError):
        SELFIE_extractor(
            assign_good_sim_files,
            'test_folder/SNR_res/'
            )
        
# tests of adj setup

def test_adj_setup_good_data(
        assign_good_sim_files
        ):
    output = SELFIE_extractor(assign_good_sim_files,
                              'test_folder/SNR_res/')
    output_2 = adj_setup(output)
    assert output_2[0] == [0.25628, 0.545628, 0.215213]
    assert output_2[-1] == [0.941004, 1.12286, 0.463877]

def test_adj_setup_missing_keyword(
        assign_good_sim_files
        ):
    output = SELFIE_extractor(assign_good_sim_files,
                              'test_folder/SNR_res/')
    del output["redshift_estimate"]
    with pytest.raises(KeyError):
        adj_setup(output)
    