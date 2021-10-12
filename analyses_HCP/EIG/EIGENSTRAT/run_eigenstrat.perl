#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}"; 
# MUST put smartpca bin directory in path for smartpca.perl to work
#
# Nhung note: must call 'module restore eig' for this script to work

$path = "/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT"; 


## train group 
$command = "smartpca.perl";
$command .= " -i $path/eigendata/hcp_cohort.geno ";
$command .= " -a $path/eigendata/hcp_cohort.snp ";
$command .= " -b $path/eigendata/hcp_split.ind ";
$command .= " -k 3 ";
$command .= " -o $path/eigen_results/train_cohort.pca ";
$command .= " -p $path/eigen_results/train_cohort.plot ";
$command .= " -e $path/eigen_results/train_cohort.eval ";
$command .= " -l $path/eigen_results/train_cohort.log ";
$command .= " -m 0 ";
$command .= " -t 3 ";
$command .= " -s 6.0 ";
$command .= " -w $path/eigendata/train_pop.txt ";
system("$command");

## test group 
$command = "smartpca.perl";
$command .= " -i $path/eigendata/hcp_cohort.geno ";
$command .= " -a $path/eigendata/hcp_cohort.snp ";
$command .= " -b $path/eigendata/hcp_split.ind ";
$command .= " -k 3 ";
$command .= " -o $path/eigen_results/test_cohort.pca ";
$command .= " -p $path/eigen_results/test_cohort.plot ";
$command .= " -e $path/eigen_results/test_cohort.eval ";
$command .= " -l $path/eigen_results/test_cohort.log ";
$command .= " -m 0 ";
$command .= " -t 3 ";
$command .= " -s 6.0 ";
$command .= " -w $path/eigendata/test_pop.txt ";
system("$command");

