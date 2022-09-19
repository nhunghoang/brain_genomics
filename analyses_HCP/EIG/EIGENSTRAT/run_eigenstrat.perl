#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}"; 
# MUST put smartpca bin directory in path for smartpca.perl to work
#
# Nhung note: must call 'module restore eig' for this script to work

$dset = "UKB"; # HCP or UKB
$path = "/data1/rubinov_lab/brain_genomics/analyses_$dset/DATA_OUTPUT"; 

## whole group 
$command = "smartpca.perl";
$command .= " -i $path/eigendata/$dset\_cohort.geno ";
$command .= " -a $path/eigendata/$dset\_cohort.snp ";
$command .= " -b $path/eigendata/$dset\_cohort.ind ";
$command .= " -k 3 ";
$command .= " -o $path/eigen_results/$dset\_cohort.pca ";
$command .= " -p $path/eigen_results/$dset\_cohort.plot ";
$command .= " -e $path/eigen_results/$dset\_cohort.eval ";
$command .= " -l $path/eigen_results/$dset\_cohort.log ";
$command .= " -m 0 ";
$command .= " -t 3 ";
$command .= " -s 6.0 ";
system("$command");

