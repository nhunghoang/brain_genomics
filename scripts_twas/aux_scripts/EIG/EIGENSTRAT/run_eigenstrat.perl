#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}"; 
# MUST put smartpca bin directory in path for smartpca.perl to work
#
# Nhung note: must call 'module restore eig' for this script to work

$path = "/data1/rubinov_lab/brain_genomics/scripts_twas/inputs_HCP"; 

## whole group 
$command = "smartpca.perl";
$command .= " -i $path/eigen\_input/HCP\_cohort.geno ";
$command .= " -a $path/eigen\_input/HCP\_cohort.snp ";
$command .= " -b $path/eigen\_input/HCP\_cohort.ind ";
$command .= " -k 3 ";
$command .= " -o $path/eigen\_output/HCP\_cohort.pca ";
$command .= " -p $path/eigen\_output/HCP\_cohort.plot ";
$command .= " -e $path/eigen\_output/HCP\_cohort.eval ";
$command .= " -l $path/eigen\_output/HCP\_cohort.log ";
$command .= " -m 0 ";
$command .= " -t 3 ";
$command .= " -s 6.0 ";
system("$command");

