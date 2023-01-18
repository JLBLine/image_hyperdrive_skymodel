# image_hyperdrive_skymodel
Scripts to directly image a hyperdrive skymodel without going into uv space. First image to a SIN projection, then optionally mosaic into a healpix projection. Final insane option to convert that healpix image back into a hyperdrive `.yaml` sky model.

## Installation
Currently relies on a specific branch (sigh) of `SHAMFI` so it can do shapelets, so get a copy of that and then change branch and `pip install`:

```
$ cd /somewhere/you/put/software
$ git clone https://github.com/JLBLine/SHAMFI
$ git checkout new_coord_system
$ pip install .
```

Then just clone this repo to your machine to get the code
```
git clone https://github.com/JLBLine/image_hyperdrive_skymodel.git
```

Hopefully the dependencies of `SHAMFI` cover this repo, if not, install some other python modules if you get an import error lol

## Create a SIN projection image

An example command looks like this, which will make a 1500 by 1500 pixel image of resolution 0.005 degrees, that gets convolved with a restoring beam of major,minor 0.05 deg. Need the high pixel resolution to correctly image the shapelet models, they need high resolution to render well. `ra_cent` and `dec_cent` obviously set the centre of the image. The catalogue contents get extrapolated to whatever `freq` you set. This will create an image called 
`test_ra000.0_dec-30.0_freq167035000.0.fits`

```
time python image_yaml_srclist.py \
    --naxis=1500 --radec_reso=0.005 \
    --bmaj=0.05 --bmin=0.05 \
    --ra_cent=0.0 --dec_cent=-30.0 \
    --bpa=0.0 \
    --srclist=srclist_pumav3_EoR0LoBESv2_fixedEoR1pietro+ForA_phase1+2_CenA-GP.yaml \
    --freq=167.035e+6 \
    --output_name=test
```

it'll print out a bunch of lines like

```
COM001613-3123 Gaussian Minor axis was tiny, setting to half pixel width
        RA,Dec 4.02,-31.40 major,minor (arcsec) 149.0,0.1, flux (Jy) 0.129
```
it makes this warning when a Gaussian component has a teeny major or minor axis, which breaks the script. Anyways congratulations, you've made an image

## Create SIN images across the entire sky
I've included a `bash` script that will make enough SIN images to cover the entire sky model that we have, for a given frequency. Example usage is like this:

```
./image_srclist.sh \
    -s srclist_pumav3_EoR0LoBESv2_fixedEoR1pietro+ForA_phase1+2_CenA-GP.yaml \
    -f 167.035e+6 \
    -o full_cat_SIN_image
```

will make a bunch of files in a folder `full_cat_SIN_image`:

```
$ ls full_cat_SIN_image/
full_cat_SIN_image_ra000.0_dec00.0_freq167035000.0.fits
full_cat_SIN_image_ra000.0_dec-30.0_freq167035000.0.fits
full_cat_SIN_image_ra000.0_dec30.0_freq167035000.0.fits
full_cat_SIN_image_ra000.0_dec-60.0_freq167035000.0.fits
full_cat_SIN_image_ra000.0_dec-90.0_freq167035000.0.fits
etc etc
```

> :warning: This makes 14GB of images and takes 2.5 hours for a single freq

## Mosaic and project to healpix
Do something like this command, which will grab everything in the `full_cat_SIN_image` dir and stick into some healpix images of the requested `nside`.

```
python project_to_healpix.py \
    --SIN_prepend=full_cat_SIN_image \
    --freq=${freq} \
    --output_prepend=fancy_mosaic \
    --nside=2048
```

This will create the files

```
fancy_mosaic_freq167035000.0_nside0256_map.fits
fancy_mosaic_freq167035000.0_nside0256_map_unweighted.fits
fancy_mosaic_freq167035000.0_nside0256_weights.fits
```

where `*map.fits` is just `*map_unweighted.fits` divided by `_weights.fits`, which is a count of all the footprints of the interpolated SIN images onto the final healpix format. This needs some thought as to what is best to do.

## Convert into yaml srclist again
I have just copied and pasted the script for the 21cm model that coverts from healpix to yaml. It's hard coded and mean. You gotta fix this Jaiden, probably need to sit down with me and ask questions. Makes a seperate model for each coarseband, which you supply as an index (zero indexed, so to make band one, you do)

```
python convert_healpix_to_yaml.py 0
```
