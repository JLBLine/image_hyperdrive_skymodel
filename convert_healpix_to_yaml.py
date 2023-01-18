import healpy as hp
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from glob import glob

def write_yaml_coarseband(healpix_cube, freqs, ra_rad, dec_rad, filename):

    flux_cut = 1e-12
    
    num_skipped = 0
    
    with open(filename, 'w') as outfile:
        outfile.write("le-sky-model:\n")
    
        for coord_ind, ra, dec in zip(np.arange(len(ra_rad)), ra_rad, dec_rad):

            ##Check if all fluxes are less than some cutoff: don't bother
            ##adding information if it's essentially a zero
            if (np.abs(healpix_cube[:, coord_ind]) < flux_cut).all():
                num_skipped += 1
            else:
                outfile.write(f"  - ra: {ra*(180/(np.pi)):.10f}\n")
                outfile.write(f"    dec: {dec*(180/(np.pi)):.10f}\n")
                outfile.write("    comp_type: point\n")

                outfile.write("    flux_type:\n")
                outfile.write("      list:\n")

                for freq_ind, freq in enumerate(freqs):
                    flux = healpix_cube[freq_ind, coord_ind]
                    outfile.write(f"        - freq: {freq}\n")
                    outfile.write(f"          i: {flux:.12f}\n")
                
    print(f"Skipped {num_skipped} out of {hp.nside2npix(2048)} ({100*(num_skipped/hp.nside2npix(2048)):.1f}%)")

    return

if __name__ == '__main__':

    from sys import argv

    nside = 2048

    coarse_band = int(argv[1])
    
    coarse_bw = 1.28e+6
    num_coarse_bands = 24
    highband_low_edge = 167.035e+6

    band_low_edges = np.arange(highband_low_edge, highband_low_edge + num_coarse_bands*coarse_bw, coarse_bw)
    
    all_files = np.array(sorted(glob("healpix_tiled/*unweighted*.fits")))
    
    all_freqs = []
    for filename in all_files:
        freq = float(filename.split('_')[-1][1:-4])
        all_freqs.append(freq)
        
    all_freqs = np.array(all_freqs)    
    edge_buffer = 80e+3
    
    band_low_edge = band_low_edges[coarse_band]
    band_high_edge = band_low_edge + coarse_bw + edge_buffer
    
    band_low_edge -= edge_buffer
    
    these_slices = np.where((band_low_edge < all_freqs) & (all_freqs < band_high_edge))
    
    num_freqs = len(these_slices[0])
    
    print(f"For band {coarse_band + 1} there are {num_freqs} slices")
    print(all_freqs[these_slices])
    
    freqs = all_freqs[these_slices]
    
    print(freqs.max(), band_high_edge)
    
    print(freqs)
    
    if freqs.max() < band_high_edge:
        freqs = np.append(freqs, freqs[-1] + (freqs[-1] - freqs[-2]))
        
    print(freqs[-2] - freqs[-1])
    
    healpix_cube = np.empty((num_freqs, hp.nside2npix(nside)))
    
    for freq_ind, filename in enumerate(all_files[these_slices]):
        this_map = hp.read_map(filename)
        print("max, min, sum", this_map.min(), this_map.max(), np.sum(this_map))
        healpix_cube[freq_ind, :] = this_map

 
        
    dec_rad, ra_rad = hp.pixelfunc.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    dec_rad = np.pi/2 - dec_rad
    
    freqs = all_freqs[these_slices]
    
    filename = f"yaml_models/21cm_skymodel_band{coarse_band+1:02d}.yaml"
    
    write_yaml_coarseband(healpix_cube, freqs, ra_rad, dec_rad, filename)
    
    ##Do one time to check you have your coord system correct
    # hp.write_map("ra_coords_healpix.fits", ra_rad/D2R, overwrite=True)
    # hp.write_map("dec_coords_healpix.fits", dec_rad/D2R, overwrite=True)
