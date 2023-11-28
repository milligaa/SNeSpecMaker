from scipy import integrate
import scipy.ndimage
import numpy as np


def point_convolute(seeing, sne_mag) -> float:
    """
    Calculates magnitude of a SNe down a 4MOST fibre by accounting for
    seeing effects.

    :param float seeing: the atmospheric seeing FWHM in arcseconds.
    :param float sne_mag: magnitude of SNe in AB-mag in LSST r-band.
    :return: effective mag of SNe down a 4MOST fibre.
    :rtype: float.
    """

    # first turn the seeing value into a sigma for the gaussian
    FWHM = seeing
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))

    # generate a radial guassian (gaussian multiplied by 2 pi x)
    gaussian = lambda x: 2 * np.pi * x * np.e ** (-(x**2) / (2 * (sigma**2)))

    # integrate to infinity to get the normalisation
    normalisation = integrate.quad(gaussian, 0, np.inf)

    # integrate again to fibre radius with normalisation
    gaussian_norm = (
        lambda x: 2
        * np.pi
        * x
        * np.e ** (-(x**2) / (2 * (sigma**2)))
        / normalisation[0]
    )
    fraction = integrate.quad(gaussian_norm, 0, 0.725)[0]

    print(fraction)

    new_mag = sne_mag - 2.5 * np.log(fraction)
    return new_mag


def effective_fibre_mag(sep, gal_ddlr, sersic_index,
                        galmag, pix_size_arcsec, seeing):

    """
    Calculates the effective magnitude of an extended host galaxy
    object seen by a 4MOST fibre placed directly on the transient
    it hosts.

    :param float sep: seperation between host and transient in arcseconds.
    :param float ddlr: directional light radius normalised by sep.
    :param float sersic_index: sersic index to be assumed for host.
    :param float galmag: AB magnitude of host in LSST r-band filter.
    :param float pix_size_arcsec: projected pixel size in arcseconds.
    :param float seeing: atmospheric seeing FWHM in arcseconds.
    :return: effective host magnitude down 4MOST fibre
    rtype: float
    """    

    dlr = gal_ddlr * sep
    gmag = galmag
    seperation = sep

    # constant here is the zero-point flux in the LSST r-band filter bandpass
    flux_gal_total = 4.69542e-6 * (np.e ** (-gmag / 2.5))

    # this flux is equal to const * pi * Ie * Re^2 (for sersic index of 0.5)
    eff_intensity = flux_gal_total / (2.8941 * np.pi * (dlr**2))

    # then calculate effective intensity from Prugneil and Simien
    bn = (2 * sersic_index) - (1 / 3) + (0.009876 * sersic_index)

    # now create the pixel array from integer number of pixels in display range
    pixel_length = pix_size_arcsec
    pixel_no = (
        int((6 * seperation) / pixel_length) + 1
    )  # plus one because always rounds down

    int_array = []
    distance_array = []
    for i in range(pixel_no):
        sub_array = []
        sub_dist = []
        for j in range(pixel_no):
            sub_array.append(0)
            sub_dist.append(
                np.sqrt((i - (pixel_no / 2)) ** 2 + (j - (pixel_no / 2)) ** 2)
                * pixel_length
            )

        int_array.append(sub_array)
        distance_array.append(sub_dist)

    # define the intensities at each pixel
    for h in range(pixel_no):
        for g in range(pixel_no):
            int_array[h][g] = np.log(
                eff_intensity
                * np.e
                ** (
                    -bn * (((distance_array[h][g]
                             / dlr
                             ) ** (1 / sersic_index)) - 1)
                    )
            )

    min_int = min(min(int_array))

    # define conv radius using seeing (FWHM -> sigma)
    conv_radius = (seeing / np.sqrt(8 * np.log(2))) / pixel_length
    new_int_array = scipy.ndimage.gaussian_filter(
        int_array, sigma=conv_radius, mode="constant", cval=min_int
    )

    pixel_sep = seperation / pixel_length
    pixel_fibre_size = 0.725 / pixel_length
    true_array_int = np.asarray(new_int_array)

    # find which pixels fall within fibre and sum intensity
    good_coords = []
    fibre_int_pix = 0
    for h in range(pixel_no):
        for g in range(pixel_no):
            pixel_centre_x = h
            pixel_centre_y = g
            pixel_distance = (pixel_centre_x -
                              ((pixel_no / 2) + pixel_sep)) ** 2 + (
                pixel_centre_y - (pixel_no / 2)
            ) ** 2

            pixel_int = true_array_int[h][g]
            if pixel_distance <= pixel_fibre_size**2:
                fibre_int_pix = fibre_int_pix + (
                    (np.e**pixel_int) * pixel_length * pixel_length
                )
                good_coords.append([h, g])
            else:
                continue

    ratio_pix = fibre_int_pix / flux_gal_total

    if ratio_pix > 1:
        ratio_pix = 1
    else:
        ratio_pix = ratio_pix

    # now calculate effective galaxy mag
    eff_mag = gmag - 2.5 * np.log10(ratio_pix)

    return eff_mag