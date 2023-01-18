from reproject import reproject_to_healpix
from astropy.io import fits
import healpy as hp
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
import argparse

coordsys = 'icrs'
nested = True
D2R = np.pi / 180.0


def do_healpix_projection(healpix_name, nside):
    with fits.open(healpix_name) as hdu:
            
            image = hdu[0].data[0,0,:,:]
            wcs = WCS(hdu[0].header).celestial
            
            ##Image is in units of Jy/beam, we want it in units of Jy/pixel
            ra_reso = hdu[0].header['CDELT1']
            dec_reso = hdu[0].header['CDELT2']
            bmaj = hdu[0].header['BMAJ']
            bmin = hdu[0].header['BMIN']
            
            ##converts to Jy/pixel I think
            beam_area_rad = (np.pi*bmaj*bmin*D2R**2) / (4*np.log(2))
            healpix_rad = hp.nside2pixarea(nside, degrees=False)

            image *= (healpix_rad / beam_area_rad)
            
            ##Do the interpolation to healpix
            healpix_data, footprint = reproject_to_healpix((image, wcs), coordsys, nside=nside, order='nearest-neighbor')
            
            ##Footprint is just weights of where the pixels have been interpolated to
            footprint[np.isnan(healpix_data)] = 0.0
            healpix_data[np.isnan(healpix_data)] = 0.0
            
    return healpix_data, footprint

def get_arg_parser():
    parser = argparse.ArgumentParser(description="I dunno")

    parser.add_argument('--nside', default=False, type=int, required=True,
        help='Nside of healpix output image')

    parser.add_argument('--SIN_prepend', default=False, required=True,
        help='The prepend of the SIN projected images to stitch into a healpix')

    parser.add_argument('--freq', default=180e+6, type=float,
        help='Frequency of the images to stitch together')

    parser.add_argument('--output_prepend', default=False, required=True,
        help='Name to prepend all the outputs with')

    return parser


if __name__ == '__main__':

    parser = get_arg_parser()
    args = parser.parse_args()

    ##Where the full map and weights get stored
    full_sky_map = np.zeros(hp.nside2npix(args.nside))
    full_sky_weights = np.zeros(hp.nside2npix(args.nside))

    ##Loop through the ra,dec images, project each one to a healpix,
    ##accumulating into one map, as well as trakcing the footprint of each
    ##image into the final map for weights
    for ra in np.arange(0., 360., 30.):
        for dec in np.arange(-60., 60., 30.):
            print(f'Projecting ra, dec, freq {ra:.1f}, {dec:.1f}, {args.freq:.6e} ')
            
            healpix_data, footprint = do_healpix_projection(f'{args.SIN_prepend}/{args.SIN_prepend}_ra{ra:05.1f}_dec{dec:04.1f}_freq{args.freq:.1f}.fits', args.nside)
            
            full_sky_map += healpix_data
            full_sky_weights += footprint
        
    ##Stick the south pole on there too
    healpix_data, footprint = do_healpix_projection(f'{args.SIN_prepend}/{args.SIN_prepend}_ra000.0_dec-90.0_freq167035000.0.fits', args.nside)

    full_sky_map += healpix_data
    full_sky_weights += footprint

    ##Do a weighted sky map - just dividing by the summed weights, nothing
    ##fancy        
    weighted_map = np.zeros(hp.nside2npix(args.nside))
    weighted_map[full_sky_weights > 0] = full_sky_map[full_sky_weights > 0] / full_sky_weights[full_sky_weights > 0]

    ##Save the raw and weighted maps, as well as the weights
    hp.write_map(f"{args.output_prepend}_freq{args.freq:.1f}_nside{args.nside:04d}_map_unweighted.fits", full_sky_map, overwrite=True)
    hp.write_map(f"{args.output_prepend}_freq{args.freq:.1f}_nside{args.nside:04d}_map.fits", weighted_map, overwrite=True)
    hp.write_map(f"{args.output_prepend}_freq{args.freq:.1f}_nside{args.nside:04d}_weights.fits", full_sky_weights, overwrite=True)
    