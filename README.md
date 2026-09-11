# SNeSpecMaker

This code should convert the pre-existing catalog-level simulations of from the
4MOST Facility Simulator into 10,000s of individually simulated spectra. These
spectra are complete with host contamination, fibre losses and the sample will
have realistically distributed properties based on the SELFIE algorithm (see
Tempel+2020a,b).

Upon pulling this repo, the could can be run using:

python SNeSpecMaker/spec_maker/generate_spectra.py -i SNeSpecMaker/spec_maker/init_input.yaml

from the directory above the SNeSpecMaker folder, although some simple edits
can fix the paths to run it from anywhere I suspect.

The inputs in the initialising .yaml file are:

spectra_save_path: Filepath pointing to where you want to save the spectra 
                   produced

host_loc: path to the host galaxy SEDs, obtained from Kinney+1996

begin: the index position in input_population to simulate spectra from

end: the index position in input_population to simulate spectra to

These two params are primarily for distributed computing or manually selecting
a small number of spectra for debugging

SNANA_SED_loc: location of SNANA-derived transient SEDs

input_population: path to transient data derived from catalog-level 4MOST
                  simulations

Both of these are created from 4MOST proprietary data. Discussion of both can
be found in Milligan+2025.