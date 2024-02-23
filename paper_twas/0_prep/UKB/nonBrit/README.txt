Preprocessing UKB subset of non-White-British individuals 
that are still in European cluster. This data is needed 
for the TWAS replication analysis. 

In actuality - the replication cohort is individuals with volume 
and genotype data who self-reported as White but are not in the 
computed White British cohort (our main 40k cohort). 

last update: Feb 2, 2024

#####################################################################################

0. Download data from server 

See brain_genomics/data_UKB/downloads/README.txt 

(At this point, we have the allPops folder with info on the 5k people who have 
vol and gen data and were not part of WhiteBrit cluster - pull self-reported 
Euro cluster from here)

#####################################################################################

1) Get eids for this cohort and their covariate data 

   0_prep/UKB/nonBrit/0_slice_cohort_eids.py 

2) Slice original bgen files down to the samples and variants of interest 

   0_prep/UKB/nonBrit/1_extract_cohort_jti.sh 
    

