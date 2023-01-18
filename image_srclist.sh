##ooh it's a help function in bash
Help()
{
   # Display Help
   echo "Makes enough SIN FITS projected images to cover the entire sky"
   echo
   echo "Syntax: image_srclist.sh -s srclist.yaml -f frequency -o output_name"
   echo "options:"
   echo "-h     Print this Help."
   echo "-s     The .yaml srclist to use as a sky model"
   echo "-f     The frequency (Hz) to extrapolate the model to"
   echo "-o     Name for the outputs and output folder"
   echo
}

##grab some options
while getopts :s:o:f:h flag
do
    case "${flag}" in
        s) srclist=${OPTARG};;
        o) output=${OPTARG};;
        f) freq=${OPTARG};;
        h) Help
           exit;;
        \?) # incorrect option
         echo "Error: Invalid option"
         Help
         exit;;
    esac
done

echo "Using srclist: ${srclist}"
echo "Extrapolating to frequency: ${freq}"
echo "Naming ouputs: ${output} "

mkdir -p ${output}

##can make theseo options if you fancy
naxis=6000
reso=0.005

start=`date +%s`

time python image_yaml_srclist.py \
            --naxis=${naxis} --radec_reso=${reso} \
            --bmaj=0.04 --bmin=0.04 \
            --ra_cent=0.0 --dec_cent=-90.0 \
            --bpa=0.0 \
            --srclist=${srclist} \
            --freq=${freq} \
            --output_name=${output}/${output}

for ra in {0..330..30}
do
    for dec in {-60..30..30}
    do
        echo "Doing ra ${ra} dec ${dec}"
        time python image_yaml_srclist.py \
            --naxis=${naxis} --radec_reso=${reso} \
            --bmaj=0.04 --bmin=0.04 \
            --ra_cent=${ra} --dec_cent=${dec} \
            --bpa=0.0 \
            --srclist=${srclist} \
            --freq=${freq} \
            --output_name=${output}/${output}
    done
done

end=`date +%s`

##Time passed in seconds
time=$((end-start))

##Time passed in minutes
time=$(awk "BEGIN {print ${time}/60.0}")

echo "Full execution time was ${time} minutes"
