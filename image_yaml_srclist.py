#!/usr/bin/env python

from astropy.convolution import convolve
from shamfi import shapelets, read_FITS_image, shapelet_coords, srclists, shamfi_plotting, image_manipulation
import numpy as np
# import matplotlib.pyplot as plt
from astropy.io import fits
# from learning_things.combine_shapelet_coefficients import gen_shape_basis_direct
from subprocess import call
# from copy import deepcopy
from os import getcwd
import argparse
from astropy.coordinates import SkyCoord

import matplotlib
# matplotlib.use('agg')

D2R = np.pi/180.0

##convert between FWHM and std dev for the gaussian function
FWHM_factor = 2. * np.sqrt(2.*np.log(2.))

def create_FITS_file(ra_cent, dec_cent, naxis, radec_reso,
                     bmaj, bmin, pa, freq, filename, data):
    hdu = fits.PrimaryHDU()

    hdu.header['NAXIS'] = 4

    hdu.header['NAXIS1'] = naxis
    hdu.header['NAXIS2'] = naxis
    hdu.header['NAXIS3'] = 1
    hdu.header['NAXIS4'] = 1

    hdu.header['CDELT1'] = -radec_reso
    hdu.header['CRPIX1'] = int(naxis // 2) + 1
    hdu.header['CRVAL1'] = ra_cent
    hdu.header['CTYPE1'] = 'RA---SIN'
    hdu.header['CUNIT1']  = 'deg'

    hdu.header['CDELT2'] = radec_reso
    hdu.header['CRPIX2'] = int(naxis // 2) + 1
    hdu.header['CRVAL2'] = dec_cent
    hdu.header['CTYPE2'] = 'DEC--SIN'
    hdu.header['CUNIT2']  = 'deg'

    hdu.header['CTYPE3'] = 'FREQ    '
    hdu.header['CRPIX3'] =  1.
    hdu.header['CRVAL3'] = freq
    hdu.header['CDELT3'] = 1e+6
    hdu.header['CUNIT3']  = 'Hz '

    hdu.header['CTYPE4'] = 'STOKES  '
    hdu.header['CRPIX4'] = 1.
    hdu.header['CRVAL4'] = 1.
    hdu.header['CDELT4'] = 1.

    hdu.header['BMAJ'] = bmaj
    hdu.header['BMIN'] = bmin
    hdu.header['BPA'] = pa
    hdu.header['BUNIT'] = 'JY/BEAM'
    hdu.header['BTYPE'] = 'Intensity'

    hdu.header['EQUINOX'] = 2000.0

    hdu.data = np.zeros((1, 1, naxis, naxis))

    hdu.data[0,0,:,:] = data

    hdu.writeto(filename, overwrite=True)

def make_all_basis_function_images(n1s, n2s, xrot, yrot, b1, b2):
    """For the given image `x_mesh, y_mesh` coords and n1,n2 pairs in `n1s,n2s`,
    create 2D image-space shapelet basis functions"""

    all_basis_functions = np.empty((len(n1s), xrot.shape[0], xrot.shape[1]))

    for ind, n1, n2 in zip(np.arange(len(n1s)),n1s, n2s):

        basis = shapelets.gen_shape_basis_direct(n1=n1,n2=n2,xrot=xrot,yrot=yrot,b1=b1,b2=b2)

        all_basis_functions[ind] = basis

    return all_basis_functions

def get_num_shape_basis(nmax):
    """For a given nmax, return the number of possible 2D shapelet basis functions"""

    return int(((nmax + 2)*(nmax + 1))/2)

def get_n1n2_index(n1, n2, nmax):
    """For a complete sequence of all n1, n2 pairs with some nmax, find the
    index of the pair in a 1D array"""

    return int(n2 + n1*(nmax -n1/2 + 1.5))

def calc_n1n2_pairs(nmax):
    """Calculates all the n1,n2 pairs used up to maximum basis order `nmax`"""

    num_basis = get_num_shape_basis(nmax)

    n1s = np.empty(num_basis)
    n2s = np.empty(num_basis)

    for n1 in range(nmax + 1):
        for n2 in range(nmax - n1 + 1):

            index = get_n1n2_index(n1, n2, nmax)

            n1s[int(index)] = n1
            n2s[int(index)] = n2

    return n1s, n2s

def make_twoD_model(coeffs, all_basis_functions):
    """`all_basis_functions` should be an array shape = (num_basis, y_len, x_len)
    containing a set of 2D basis function results over a grid of
    shape = (y_len, x_len). `coeffs` should be coefficients corresponding to
    all basis functions. Returns the 2D summed model of the coefficients
    multiplied by the basis functions"""

    mult_by_coeffs = all_basis_functions*coeffs[:, np.newaxis, np.newaxis]

    return np.sum(mult_by_coeffs, axis=0)


def do_point_source(fits_data, component, total_summed_model):
    """Add a basic point source model; just adds all the flux to a single
    pixel"""
    x,y = fits_data.wcs.all_world2pix(component.ra/D2R, component.dec/D2R, 0)
    x, y = int(np.round(x)), int(np.round(y))

    if x < 0 or x >= fits_data.len1 or y < 0 or y >= fits_data.len2:
        pass
    else:
        flux = srclists.extrapolate_component_flux(component,  args.freq)
        total_summed_model[y,x] += flux / fits_data.convert2pixel
    return total_summed_model

def main(args):

    ##Make a fake FITS file to read in via the SHAMFI machinery=================

    some_data = np.arange(args.naxis*args.naxis)
    some_data.shape = (args.naxis, args.naxis)

    freq = args.freq

    create_FITS_file(args.ra_cent, args.dec_cent, args.naxis, args.radec_reso,
                           args.bmaj, args.bmin, args.bpa, freq, 'temp_shamfi_model.fits',
                           some_data)

    fits_data = read_FITS_image.FITSInformation(f"{getcwd()}/temp_shamfi_model.fits")
    fits_data.get_radec_edgepad(edge_pad=False)
    fits_data.create_restoring_kernel()

    ##Read in the fake FITS file and setup to shapelet coord system=============

    shpcoords = shapelet_coords.ShapeletCoords(fits_data)
    shpcoords.find_good_pixels()

    # components = read_woden_srclist(args.srclist)
    components = srclists.read_hyperdrive_srclist(args.srclist)

    total_summed_model = np.zeros((args.naxis,args.naxis))

    all_ras = [component.ra for component in components]
    all_decs = [component.dec for component in components]

    coords = SkyCoord(all_ras, all_decs, unit='rad')

    inside_image = np.where(coords.contained_by(fits_data.wcs) == True)[0]

    print(f"There are {len(inside_image)} components inside the image to render")

    ##coords of the image
    x_pixel_mesh, y_pixel_mesh = np.meshgrid(np.arange(fits_data.len1),
                                             np.arange(fits_data.len2))

    ras_mesh, decs_mesh = fits_data.wcs.all_pix2world(x_pixel_mesh, y_pixel_mesh, 0)

    ras_mesh *= D2R
    decs_mesh *= D2R

    for comp_ind, component in enumerate(components[inside_image]):

        if component.comp_type == 'SHAPELET':
            # print("Doing Shapelet")

            x,y = fits_data.wcs.all_world2pix(component.ra/D2R, component.dec/D2R, 0)

            x,y = int(x), int(y)

            big_axis = np.max([component.major, component.minor])/D2R

            num_axis = 10
            half_width_box = (num_axis*big_axis) / args.radec_reso

            x_low = int(np.floor(x - half_width_box))
            x_high = int(np.ceil(x + half_width_box))
            y_low = int(np.floor(y - half_width_box))
            y_high = int(np.ceil(y + half_width_box))

            if x_low < 0: x_low = 0
            if x_high >= fits_data.len1: x_high = fits_data.len1 - 1

            if y_low < 0: y_low = 0
            if y_high >= fits_data.len2: y_high = fits_data.len2 - 1

            x_width = x_high - x_low + 1
            y_width = y_high - y_low + 1

            # print(half_width_box)

            fit_box = f"{x_low},{x_high},{y_low},{y_high}"

            # print(fit_box)

            shpcoords = shapelet_coords.ShapeletCoords(fits_data)
            shpcoords.find_good_pixels(fit_box=fit_box)
            #
            shpcoords._set_given_radec_to_zero_pixel(component.ra, component.dec)
            shpcoords.pa = component.pa

            xrot, yrot = shpcoords.radec2xy(component.major, component.minor, crop=True)

            xrot.shape = (y_width, x_width)
            yrot.shape = (y_width, x_width)

            all_basis_functions = make_all_basis_function_images(component.n1s,
                                                                 component.n2s,
                                                                 xrot, yrot,
                                                                 component.major,
                                                                 component.minor)

            ##Make the output model Jy/beam. The stored coeffs inside the srclist
            ##are normalised, so need to multiply by `flux`
            ##Need to account for pixel area to get into Jy, and convert from per pixel
            ##into per beam
            scale = fits_data.pix_area_rad / fits_data.convert2pixel
            # flux = component.calc_flux(args.freq)
            flux = srclists.extrapolate_component_flux(component,  args.freq)

            summed_model = scale*flux*make_twoD_model(component.shapelet_coeffs,
                                                      all_basis_functions)

            total_summed_model[y_low:y_high+1, x_low:x_high+1] += summed_model

        if component.comp_type == 'GAUSSIAN':
            
            small_axes = args.radec_reso*D2R / 4
            conv_2_arcsec = 3600 / D2R
            
            ##If both the major and minor axes are less than the resolution of
            ##the image this is basically a point source, so just add as a point source
            if component.major < small_axes and component.minor < small_axes:
                print(f"{component.source_name} both axes were tiny, converting Gaussian to point source")
                print(f"RA,Dec {component.ra:.2f}/D2R,{component.dec:.2f}/D2R major,minor (arcsec) {component.major*conv_2_arcsec:.1f},{component.minor*conv_2_arcsec:.1f}")
                
                total_summed_model = do_point_source(fits_data, component, total_summed_model)
            
            else:
                ##If one of the major or minor axis is tiny, we need to switch
                ##to the pixel resolution otherwise the gridding goes fucky
                
                if component.major < small_axes / 2:
                    major = small_axes / 2
                    print(f"{component.source_name} Gaussian Major axis was tiny, setting to half pixel width")
                    print(f"RA,Dec {component.ra:.2f}/D2R,{component.dec:.2f}/D2R major,minor (arcsec) {component.major*conv_2_arcsec:.1f},{component.minor*conv_2_arcsec:.1f}")
                else:
                    major = component.major
                    
                    
                if component.minor < small_axes / 2:
                    minor = small_axes / 2
                    print(f"{component.source_name} Gaussian Minor axis was tiny, setting to half pixel width")
                    print(f"\tRA,Dec {component.ra/D2R:.2f},{component.dec/D2R:.2f} major,minor (arcsec) {component.major*conv_2_arcsec:.1f},{component.minor*conv_2_arcsec:.1f}, flux (Jy) {srclists.extrapolate_component_flux(component,  args.freq):.3f}")
                else:
                    minor = component.minor
                    
                # major = component.major
                # minor = component.minor

                ra_stddev = major / (FWHM_factor)
                dec_stddev = minor / (FWHM_factor) 

                x,y = fits_data.wcs.all_world2pix(component.ra/D2R, component.dec/D2R, 0)

                x,y = int(x), int(y)

                big_axis = np.max([major, minor])/D2R

                num_axis = 5
                half_width_box = (num_axis*big_axis) / args.radec_reso

                x_low = int(np.floor(x - half_width_box))
                x_high = int(np.ceil(x + half_width_box))
                y_low = int(np.floor(y - half_width_box))
                y_high = int(np.ceil(y + half_width_box))

                if x_low < 0: x_low = 0
                if x_high >= fits_data.len1: x_high = fits_data.len1 - 1

                if y_low < 0: y_low = 0
                if y_high >= fits_data.len2: y_high = fits_data.len2 - 1

                x_width = x_high - x_low + 1
                y_width = y_high - y_low + 1

                # print(half_width_box)

                fit_box = f"{x_low},{x_high},{y_low},{y_high}"

                shpcoords = shapelet_coords.ShapeletCoords(fits_data)
                shpcoords.find_good_pixels(fit_box=fit_box)
                shpcoords._set_given_radec_to_zero_pixel(component.ra, component.dec)
                shpcoords.pa = component.pa

                xrot, yrot = shpcoords.radec2xy(major, minor, crop=True)

                xrot.shape = (y_width, x_width)
                yrot.shape = (y_width, x_width)
                
                flux = srclists.extrapolate_component_flux(component,  args.freq)

                ##We've scaled the coord system to the gaussian maj/min, so
                ##make basic gaussian function and shove in the 
                gauss_func = image_manipulation.Gaussian2D(amplitude=flux,
                                        x_mean=0.0, y_mean=0.0,
                                        x_stddev=1, y_stddev=1,
                                        theta=0.0)

                gauss_model = gauss_func(xrot, yrot)
                
                # print(gauss_model)

                if gauss_model.sum() == 0.0:
                    print(component.source_name)
                    print("Something went wrong with gaussian", component.ra, component.dec, flux)
                    print(component.major, component.minor)

                else:

                    ##Get the convertsion from Jy/beam to Jy/pixel
                    ##Scale the gaussian to subtract to match the desired integrated flux
                    # print(gauss_model.sum())
                    gauss_model *= flux / (gauss_model.sum()*fits_data.convert2pixel)

                    total_summed_model[y_low:y_high+1, x_low:x_high+1] += gauss_model

        if component.comp_type == 'POINT':
            total_summed_model = do_point_source(fits_data, component, total_summed_model)
            
        else:
            pass
            


    ##Convolve by the restoring beam to match Jy/beam
    convolved_model = convolve(total_summed_model, fits_data.rest_gauss_kern)

    create_FITS_file(args.ra_cent, args.dec_cent, args.naxis, args.radec_reso,
                     args.bmaj, args.bmin, args.bpa, freq, f'{args.output_name}_ra{args.ra_cent:05.1f}_dec{args.dec_cent:04.1f}_freq{args.freq}.fits',
                      convolved_model)

    call(f"rm {getcwd()}/temp_shamfi_model.fits", shell=True)


def get_arg_parser():
    parser = argparse.ArgumentParser(description="I dunno")

    parser.add_argument('--naxis', default=False, type=int,
        help='naxis of output image')

    parser.add_argument('--radec_reso', default=False, type=float,
        help='RA/Dec resolution of output image')

    parser.add_argument('--ra_cent', default=False, type=float,
        help='RA (deg) at image centre')

    parser.add_argument('--dec_cent', default=False, type=float,
        help='Dec (deg) at image centre')

    parser.add_argument('--bmaj', default=False, type=float,
        help='Major axis restoring beam')

    parser.add_argument('--bmin', default=False, type=float,
        help='Minor axis restoring beam')

    parser.add_argument('--bpa', default=False, type=float,
        help='PA restoring beam')

    parser.add_argument('--freq', default=180e+6, type=float,
        help='Frequency to scale the image to')

    parser.add_argument('--srclist', default=False,
        help='path to srclist to image')

    parser.add_argument('--output_name', default=False,
        help='Name for the outputs')

    return parser


if __name__ == '__main__':

    parser = get_arg_parser()
    args = parser.parse_args()

    main(args)
