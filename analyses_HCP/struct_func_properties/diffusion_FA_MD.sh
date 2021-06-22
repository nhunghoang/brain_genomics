#!/bin/bash

## Script for calculating structural-diffusion properties like: fractional anisotropy (FA) and mean diffusivity (MD) from Diffusion-weighted MRI data
##this bash file use FSL dtifit command through HCP subjects directory 
## get raw diffusion data + bvecs + bvals + brain mask, to generate FA and MD as outputs  

## Neda SP, June 2021 


module load GCC/5.4.0-2.26 OpenMPI/1.10.3 FSL/5.0.10

data_dir="/data1/rubinov_lab/Neda/data" 

for dir in ${data_dir}/*/       ## list directories in the form "/data1/rubinov_lab/Neda/data/subjids/"
do
[[ "${dir}" == */ ]] && dir="${dir: : -1}"   ## remove the trailing "/"

cd $dir/T1w/Diffusion \
echo "$(pwd)"
## echo "$dir/T1w/Diffusion/data.nii.gz"
 dtifit -k data.nii.gz -r bvecs -b bvals -m nodif_brain_mask.nii.gz -o dti_ \

echo "Done!"

##output: dti__L1.nii.gz  dti__MO.nii.gz  dti__V3.nii.gz
##        dti__L2.nii.gz  dti__S0.nii.gz  dti__L3.nii.gz 
##        dti__V1.nii.gz  dti__FA.nii.gz  dti__MD.nii.gz  dti__V2.nii.gz 

done
