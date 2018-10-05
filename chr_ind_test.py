# def amr_snps(vcf_file, chr, thresh):

import os
import allel
import os.path

import h5py
import numpy as np

# file source and destination for i/o
SOURCE_PATH = 'C:'+os.sep+'Users'+os.sep+'rsoto'+os.sep+'Downloads'+os.sep
DEST_PATH = 'D:'+os.sep+'UPRM'+os.sep+'Genomes'+os.sep

# vcf source file name
VCF_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

# converted file name
HDF_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.h5'
CSV_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.csv'

# if the hdf5 file does not exist yet, start the conversion
if not os.path.exists(DEST_PATH+HDF_FILE):
    allel.vcf_to_hdf5(SOURCE_PATH+VCF_FILE, DEST_PATH+HDF_FILE, fields='*', overwrite=True)

# load the data from the hdf5 file
chr22 = h5py.File(DEST_PATH+HDF_FILE, 'r')

snp_id = chr22['variants/ID']
gt = chr22['calldata/GT']
all_inds = chr22['samples']

# load allele frequencies for each population
afr_af = np.array(chr22['variants/AFR_AF'], dtype=np.float64)
amr_af = np.array(chr22['variants/AMR_AF'], dtype=np.float64)
sas_af = np.array(chr22['variants/SAS_AF'], dtype=np.float64)
eas_af = np.array(chr22['variants/EAS_AF'], dtype=np.float64)
eur_af = np.array(chr22['variants/EUR_AF'], dtype=np.float64)

# create empty arrays
sum_af = np.zeros([np.size(afr_af[:,1]),1], dtype=np.float64)   # sum of AF
af_ind = np.zeros([np.size(afr_af[:,1]),1], dtype=np.float64)      # retain original indices

# concat population arrays, excluding AMR
no_amr_af = np.concatenate((af_ind, afr_af, sas_af, eas_af, eur_af), axis=1)

# set maximum sum of AF threshold
# if the sum exceeds thresh, we will not consider that SNP
thresh = 0

# cycle through all SNPs, appending the SNP
amr_snps = []
for i, af in enumerate(no_amr_af):
    if np.all(af == np.isnan):
        continue
    elif np.nansum(af) > thresh:
        no_amr_af[i,:] = np.nan
    elif np.nansum(af) <= thresh:
        amr_snps.append(i)

# all_gts = allel.GenotypeArray(gt[amr_snps, 0:1])  # get genotypes for AMR SNPs
# allel.GenotypeArray(gt[:,0:1])                    # count: 1103547

# read PUR individual list from tsv file (found on 1000Genomes)
PUR_ids = np.genfromtxt(DEST_PATH+os.sep+"PUR_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]

# set up lists of individuals for reference
nonPUR_gt = []                 # nonPUR
PUR_gt = []                    # PUR
for i, ind in enumerate(all_inds):
    for j, pur_ind in enumerate(PUR_ids):
        if ind == pur_ind:
            PUR_gt.append(i)
        else:
            if i not in ALL_gt:
                nonPUR_gt.append(i)
            else:
                continue

print("\n\nThere are "+str(len(PUR_gt))+" PUR individuals in this dataset")

# pull genotypes of Puerto Rican individuals (HGXXXXXX)
PUR_inds = allel.GenotypeArray(gt[:,PUR_gt])
nonPUR_inds = allel.GenotypeArray(gt[:,nonPUR_gt])

# eliminate non-AMR genotypes in Puerto Rican individuals
# PUR_AMR_gts = allel.GenotypeArray(PUR_inds[amr_snps,:])

# iterate over all PUR individuals to eliminate variants not found in PUR
# WRONG: How do we determine this?
# Do we go by frequency in a population vs world population?
# If so, I have to rethink all of this
# freq = (occurences of a combination)/(# of samples)

PUR_snps = []
for i, ind in enumerate(PUR_inds):
    if np.all(PUR_inds[i,:,0] == 0) and np.all(PUR_inds[i,:,1] == 0):
        continue
    else:
        PUR_snps.append(i)

# MEETING COMMENTS (OCT 4):
# # 0 is usually more frequent than 1
# # if alternate allele is not present
# # Take files and eliminate all SNPs with more than 1 alternate allele
# # For cases where most people are 1|1
# # Ideas:
# # # In PUR: sum > 0; in world = 0
# # # 1) If any genotype has allele > 1; eliminate it
# # # 2) Eliminate all SNPS where not (see up)
# # See distribution of all SNPs along the chromosome
# # One SNP every 100,000 bases
# # We expect the distance between one SNP and another to be 100,000 bases
# CELSR1 gene ancestry probably wrong
#
# NEXT:
# # Identify distance between SNPs
# # Reverse deconvolution of ancestry data
# # Do our own deconvolution to see if it matches 1000Genomes' deconvolution
# # # If our model coincides with theirs, confidence is improved
# # # Should be different from 1000Genomes'
# Deconvolution Ideas
# # Eliminate ASW (they're African American; might bias our results)
# # # ASW was analyzed to calculate AFR AF
# # Isolate all EUR and AFR SNPs (similar to PUR SNP isolation)
# # # We expect 1 SNP for every 10,000 bases
# # # If there are 3 EUR SNPs where there are PUR SNPs, that area might be of EUR origin
# # # # We should look at the area around the PUR SNPS
# # # Based only on the unique SNPs for AFR and EUR


print("\nThere are "+str(len(PUR_snps))+" SNPs present in PUR, out of "+str(len(af_ind)))

chr22.close()
