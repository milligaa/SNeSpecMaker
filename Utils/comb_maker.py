import astropy.units as u
from astropy.io import ascii, fits
from qmostetc import QMostObservatory, Spectrum, L1DXU
from scipy import interpolate
from astropy.table import Table
from sigfig import round
import numpy as np

def Comb_Maker(SNe_data: str, Gal_data: str, Gal_mag: float, SNe_mag: float,
               Redshift: float, Gal_type: str, SN_type: str, texp: float,
               seeing: float, spec_save_path: str):
    """
    This code takes all of the parameters from the Selfie and population simulations and uses then to scale and
    combine SNe and Host templates. It also calculates the combined Host/SNe SNR.
    :param str SNe_data: filepath leading to the SNe template to be used.
    :param str Gal_data: filepath leading to the Galaxy template to be used.
    :param float Gal_mag: galaxy magnitude to be used in AB mag in LSST r-band.
    :param float SNe_mag: SNe magnitude to be used in AB mag in LSST r-band.
    :param float Redshift: redshift of contaminated spectrum.
    :param str Gal_type: string indicating Host species.
    :param str SN_type: string indicating SNe species.
    :param float texp: exposure time for spectrum, for ETC exposure in minutes.
    :param float
    """

    #defining parameters for use down the code
    G_mag = Gal_mag * u.ABmag
    S_mag = SNe_mag * u.ABmag
    redshift = Redshift


    #importing the data to be used
    spectrum = fits.open(SNe_data)
    spec_data = spectrum[1].data

    #need to correct these to redshift = 0 from the redshift listed in the filenames
    str_begin = SNe_data.find('shift') + 5
    str_end = SNe_data.find('.fits')
    redshoft = float(SNe_data[str_begin:str_end])

    lam = []
    intens = []
    for i in range(len(spec_data)):
        lam.append(spec_data[i][0] / (1 + redshoft))
        intens.append(spec_data[i][1])  

    #making a pectrum object and then using the magnitude of the sectrum to set a new magnitude

    wavelength = lam * u.nm / 10
    photon_flux = intens * u.erg / (u.cm ** 2 * u.s * u.angstrom)
    wavelength_z = wavelength * (1 + redshift)

    spec = Spectrum(wavelength_z, photon_flux)
    spec_mag = spec.get_mag(u.ABmag, 'LSST_LSST.r')

    #note the 0.545 mag factor added to Smag to account for fibre loss
    flux_ratio = 10 ** ((S_mag.value - spec_mag.value)/(-2.5)) 

    new_flux = photon_flux * flux_ratio

    #okie dokie now lets try and repeat the above for the template data

    f2=open(Gal_data)
    lines = f2.readlines()
    lam2 = []
    intens2 = [] 
    for line in lines:
        p = line.split()
        lam2.append(float(p[0]))
        intens2.append(float(p[1]))
    f2.close()

    #only take a part of the spectrum, too wide a range breaks the interpolation later
    wavelength2 = lam2 * u.nm / 10     #divide by ten to 'make' wavelength nm and agree with units of SNe data
    photon_flux2 = intens2 * u.erg / (u.cm ** 2 * u.s * u.angstrom)
    wavelength_z2 = wavelength2 * (1 + redshift)

    spec2 = Spectrum(wavelength_z2, photon_flux2)
    spec_mag2 = spec2.get_mag(u.ABmag, 'LSST_LSST.r')

    flux_ratio2 = 10 ** ((G_mag.value - spec_mag2.value)/(-2.5)) 

    new_flux2 = photon_flux2 * flux_ratio2

    #now we interpolate the galaxy spectrum on to the same wavelength values as the SNe spectrum
    #we start by defining the redshifted spectrum eavelengths here

    g_wave_z = spec2.wavelength
    g_flux_value = new_flux2.value
    g_wave_value = g_wave_z.value        

    #now we perform the interpolation

    interp = interpolate.splrep(g_wave_value, g_flux_value)
    interp_gal = interpolate.splev(wavelength_z.value, interp)


    #now we add varing amounts of contamination from galaxy to SN
    contaminated_flux = []

    for i in range(len(lam)):
        interp_gal[i] = interp_gal[i]
        contaminated_flux.append(new_flux[i].value + interp_gal[i])


    
    # Object to simulate the 4MOST observatory, including atmosphere,
    # telescope, spectrograph, CCD.
    qmost = QMostObservatory('lrs')
    obs = qmost(33.55730976*u.deg, seeing*u.arcsec, 'dark')

    #define a new spectrum obk=ject with the redshifted wavelength and the SN flux with host galaxy contamination and set up correct spectral units
    new_spec = Spectrum(wavelength_z, contaminated_flux)
    new_spec.flux= new_spec.flux * u.erg / (u.cm ** 2 * u.s * u.angstrom)
    new_mag = new_spec.get_mag(u.ABmag, 'LSST_LSST.r')
    #divide by arcseconds squared for flat spatial distro, already account for fibre size by adding 0.545 mag to SNe
    new_spec.flux = new_spec.flux / (1.665 * u.arcsec * u.arcsec)

    #set the new spectrum as the target for observation and then set the observational parameter

    obs.set_target(new_spec, 'flat')
    res = obs.expose(texp*u.min, 1) 

    #hopefully this little bit of code here will produce the out put but with all 3 arms joined
    #we create a wrapper for the output (called dxu here) and then write the joined data to a fits file

    dxu = L1DXU(qmost, res, texp*u.min)   #doesn't like exposures sperarated into blocks, not sure if matters
    comb_spec = dxu.joined_spectrum()

        
    #first import the transmission spectrum
    T_data = '/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/adjusted_flat_spec.txt'
    T_data_table = Table.read(T_data, format = 'csv', delimiter = ' ')

    lam_t = []
    intens_t = []
    for rty in range(len(T_data_table)):
        lam_t.append(T_data_table[rty][0])
        intens_t.append(T_data_table[rty][1])

    #for now the T spectrum is generated using a combined spectrum so wavelength array matches without need for interpolation
    #THIS MAY NOT ALWAYS BE TRUE, BE WARY TRAVELLER
    corrected_flux = []
    for q in range(len(lam_t)):
        corrected_flux.append(comb_spec['FLUX'][q].value / intens_t[q])

    comb_spec['FLUX'][:] = corrected_flux * u.erg / (u.cm ** 2 * u.s * u.angstrom)

    #setting up strings for data wanted in the file names for writing the combined spec
    gt = str(Gal_type)
    st = str(SN_type)
    sm = str(round(SNe_mag, 5))
    gm = str(round(Gal_mag, 5))
    re = str(round(Redshift, 5))
    exp = str(round(texp, 5))
    
    #exreacting the data from the combined spectra

    comb_table = Table()
    comb_table['wavelength'] = comb_spec['WAVE']
    comb_table['flux'] = comb_spec['FLUX']
    comb_table['err'] = comb_spec['ERR_FLUX']
    comb_table['flux_nos'] = comb_spec['FLUX_NOSS']
    comb_table['err_nos'] = comb_spec['ERR_FLUX_NOSS']

    #setting up the system path for saving the file in the right place and then writing it there
    
    name_of_file = spec_save_path+gt+st+'Smag'+sm+'Gmag'+gm+'z'+re+'texp'+exp+'.txt'
    
    ascii.write(comb_table, name_of_file, format = 'no_header', overwrite = False)
    
    #calculating the combined galaxy and SNe SNR for the contaminated spectrum
    comb_SNR_unbin = []
    for no in range(len(comb_spec['WAVE'])):
        comb_SNR_unbin.append(comb_spec['FLUX'][no]/comb_spec['ERR_FLUX_NOSS'][no]) 


    #print(comb_spec['WAVE'])
    snr_bin = []
    noise_bin = []
    target_bin = []
    wl_bin = []
    bin_number = int((17315 - 3315)/ 60)
    #print('this is bin number', bin_number)
    for ijk in range(bin_number):
        noise_sum = 0
        target_sum = 0
        wl_sum = 0
        tot = 0
        for jki in range(60):
            target_sum = target_sum + comb_spec['FLUX'][3315 + (ijk * 60) + jki].value
            noise_sum = noise_sum + (comb_spec['ERR_FLUX_NOSS'][3315 + (ijk * 60) + jki].value ** 2)
            wl_sum = wl_sum + comb_spec['WAVE'][3315 + (ijk * 60) + jki].value
            tot = tot + 1

        snr_bin.append(target_sum / np.sqrt(noise_sum))
        noise_bin.append(np.sqrt(noise_sum) / tot)
        target_bin.append(target_sum / tot)
        wl_bin.append(wl_sum / tot)

    comb_snr = np.sum(snr_bin) / bin_number
    print('binned TiDES SNR from L1 output:' , comb_snr)

    return [comb_snr, new_mag, comb_spec['WAVE'], comb_spec['FLUX']]