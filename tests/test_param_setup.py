import pytest
from astropy.table import Table
from SNeSpecMaker.Utils.param_setup import assign_sim_files, SELFIE_extractor, adj_setup

good_data_path = 'SNeSpecMaker/tests/param_setup_test_data/good_test_data/'
missing_file_path = 'SNeSpecMaker/tests/param_setup_test_data/test_data_missing_file/'
wrong_filetype_path = 'SNeSpecMaker/tests/param_setup_test_data/wrong_filetype.txt'
csv_wrong_header_path = 'SNeSpecMaker/tests/param_setup_test_data/csv_wrong_header.csv'
fits_wrong_header_path = 'SNeSpecMaker/tests/param_setup_test_data/fits_wrong_header.fits'
non_text_data_path = 'SNeSpecMaker/tests/param_setup_test_data/test_non_text_data.png'
S238_S1001_missing_col = 'SNeSpecMaker/tests/param_setup_test_data/S238_S1001_test_missing_col.csv'
dump_path = 'SNeSpecMaker/tests/param_setup_test_data/dump/'
adj_setup_input = 'SNeSpecMaker/tests/param_setup_test_data/good_test_data/good_pop_data.csv'

#this filepath is required to test adj_setup, edit this as needed
SNANA_temp_path = '/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/SNANA_temps/'

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
        dump_path
        )
    
    assert len(table.columns) == 31
    assert len(table) == 3
    assert table[1]['redshift_estimate'] == 0.545628


def test_extract_empty_dictionary_csv(assign_good_sim_files):
    
    assign_good_sim_files["galpop"] = ''

    with pytest.raises(FileNotFoundError):
        SELFIE_extractor(
            assign_good_sim_files,
            dump_path
            )


def test_extract_empty_dictionary_fits(assign_good_sim_files):
    
    assign_good_sim_files["fibres"] = ''

    with pytest.raises(FileNotFoundError):
        SELFIE_extractor(
            assign_good_sim_files,
            dump_path
            )

def test_extract_csv_wrong_columns(assign_good_sim_files):

    assign_good_sim_files["selfie"] = S238_S1001_missing_col

    with pytest.raises(KeyError):
        SELFIE_extractor(
            assign_good_sim_files,
            dump_path
            )
        
# tests of adj setup

@pytest.fixture
def input_adj_setup():
    data_to_use = Table.read(adj_setup_input,
                                 format = 'csv', delimiter = ',')

    return data_to_use

def test_adj_setup_good_data(input_adj_setup):

    output = adj_setup(input_adj_setup, SNANA_temp_path)
    print(output[0])
    assert output[9] == [17535717, 22308419, 27080102]
    assert output[-1] == ['II', 'SALT2.WF', 'II']

def test_adj_setup_missing_keyword(
        input_adj_setup
        ):
    input = input_adj_setup
    del input["redshift_estimate"]
    with pytest.raises(KeyError):
        adj_setup(input, SNANA_temp_path)

#test of assign_host

#test of get_SN_and_host_type