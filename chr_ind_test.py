def plot_snps(PUR_snps, dest):

    import h5py
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    idx = [i[0] for i in PUR_snps]
    snp_num = [i[1] for i in PUR_snps]

    sorted_PUR_snps = sorted(PUR_snps, key = lambda x: pos[x[0]])

    total_dist = max(pos[idx]) - min(pos[idx])

    spatial_binnr = int(total_dist/1686)

    snp_binnr = 1103547//100

    print()
    print("Plotting putative Puerto Rican SNPs...")

    plt.figure()
    count_dist = sns.distplot(pos[idx], bins=spatial_binnr, kde=False, norm_hist=False)
    count_dist.set_title("Raw counts distribution")
    plt.savefig(dest+'count_dist_chr22.png')

    plt.figure()
    norm_dist = sns.distplot(pos[idx], bins=spatial_binnr, norm_hist=True)
    norm_dist.set_title("Normalized distribution")
    plt.savefig(dest+'norm_dist_chr22.png')

    plt.figure()
    all_snps_count = sns.distplot(np.where(snp_id_cp[1,:] == 1), bins=spatial_binnr, kde=False, norm_hist=False)
    all_snps_count.set_title("Raw counts distribution of all SNPs")
    plt.savefig(dest+'raw_counts_snp_dist_chr22.png')

# def amr_snps(vcf_file, chrom, thresh):

import os
import io
import ast
import allel
import os.path
import itertools

import h5py
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import time        # to time the program
import datetime    # to print time at start of program
print("Time at start:\t"+str(datetime.datetime.now().time()))

start = time.time()

# file source and destination for i/o
SOURCE_PATH = 'C:'+os.sep+'Users'+os.sep+'rsoto'+os.sep+'Downloads'+os.sep
DEST_PATH = 'D:'+os.sep+'UPRM'+os.sep+'Genomes'+os.sep

# vcf source file name
VCF_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

# converted file name
HDF_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.h5'
# CSV_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.csv'

# if the hdf5 file does not exist yet, start the conversion
if not os.path.exists(DEST_PATH+HDF_FILE):
    print("\nHDF file does not exist. Creating one...")
    allel.vcf_to_hdf5(SOURCE_PATH+VCF_FILE, DEST_PATH+HDF_FILE, fields='*', overwrite=True)

# load the data from the hdf5 file
print("\nLoading HDF file...")
chr22 = h5py.File(DEST_PATH+HDF_FILE, 'r')
chrom = 22

snp_id = chr22['variants/ID'].value
z_arr = np.zeros((1, np.size(snp_id)))
snp_id_cp = np.vstack((snp_id, z_arr))     # add row full of zeros for future bool masking

pos = chr22['variants/POS'].value
gt = chr22['calldata/GT']
all_inds = chr22['samples']

# load allele frequencies for each population
# afr_af = np.array(chr22['variants/AFR_AF'], dtype=np.float64)
# amr_af = np.array(chr22['variants/AMR_AF'], dtype=np.float64)
# sas_af = np.array(chr22['variants/SAS_AF'], dtype=np.float64)
# eas_af = np.array(chr22['variants/EAS_AF'], dtype=np.float64)
# eur_af = np.array(chr22['variants/EUR_AF'], dtype=np.float64)

# create empty arrays
# sum_af = np.zeros([np.size(afr_af[:,1]),1], dtype=np.float64)   # sum of AF
# af_ind = np.zeros([np.size(afr_af[:,1]),1], dtype=np.float64)      # retain original indices

# concat population arrays, excluding AMR
# no_amr_af = np.concatenate((af_ind, afr_af, sas_af, eas_af, eur_af), axis=1)

# set maximum sum of AF threshold
# if the sum exceeds thresh, we will not consider that SNP
# thresh = 0                 # count: 1103547

# read PUR individual list from tsv file (found on 1000Genomes)
PUR_ids = np.genfromtxt(DEST_PATH+os.sep+"PUR_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]

if not os.path.exists(DEST_PATH+'chr'+str(chrom)+'_PUR_inds.txt'):
   with open(DEST_PATH+'chr'+str(chrom)+'_PUR_inds.txt', 'w') as ff:
       ff.write(str(PUR_ids))

# set up lists of individuals for reference
nonPUR_gt = []                 # nonPUR
PUR_gt = []                    # PUR
for i, ind in enumerate(all_inds):
    for j, pur_ind in enumerate(PUR_ids):
        if ind == pur_ind:
            PUR_gt.append(i)

for i, ind in enumerate(all_inds):
    if ind not in PUR_ids:
        nonPUR_gt.append(i)

# Check whether calculations were correct
print("There are "+str(len(PUR_gt))+" PUR individuals in this dataset")
print("There are "+str(len(nonPUR_gt))+" non-PUR individuals in this dataset")

# pull genotypes of Puerto Rican individuals (HGXXXXXX)
# 3D array: [variants, individuals, genotype]
print()
print("Creating genotype arrays...")
PUR_inds = allel.GenotypeArray(gt[:,PUR_gt])
nonPUR_inds = allel.GenotypeArray(gt[:,nonPUR_gt])

#if not os.path.exists(DEST_PATH+'chr'+str(chrom)+'PUR_SNPs.txt'):
print()
print("Creating SNP file")
print("Analyzing SNPs:")
tmpPUR_snps = []
nonPUR_snps = []
count = 0
for i, snp in enumerate(snp_id):
    if count == 10000:
        print("Analyzing SNP " + str(i))
        count = 0
    count += 1
    if np.any(np.nansum(PUR_inds[i,:,:], axis = 0)) > 0 and np.any(np.nansum(nonPUR_inds[i,:,:], axis = 0)) == 0:
        tmpPUR_snps.append(tuple([i,snp]))
        snp_id_cp[1,i] = 1
    elif np.any(np.nansum(PUR_inds[i,:,:], axis = 0)) < 208 and np.any(np.nansum(nonPUR_inds[i,:,:], axis = 0)) == 2400:
        tmpPUR_snps.append(tuple([i,snp]))
        snp_id_cp[1,i] = 1
    else:
        nonPUR_snps.append(tuple([i,snp]))

tmpPUR_snps = set(tmpPUR_snps)
nonPUR_snps = set(nonPUR_snps)

print()
print("Filtering SNPs")
count = 0
PUR_snps = []
for j,prsnp in enumerate(tmpPUR_snps):
    if count == 1000:
        print("Filtering SNP "+str(j))
        count = 0
    count += 1
    if prsnp not in nonPUR_snps:
        PUR_snps.append(prsnp)

#    with open(DEST_PATH+'chr'+str(chrom)+'PUR_SNPs.txt', 'w') as ff:
#        ff.write(str(PUR_snps))

#else:
#    print()
#    print("Loading PUR SNPs")
#    PUR_snps = []
#    with open(DEST_PATH+'chr'+str(chrom)+'PUR_SNPs.txt', 'r') as ff:
#        PUR_snps = ast.literal_eval(ff.read())

print()
print("PUR SNPs found:\t"+str(len(PUR_snps)))
print("Which is "+str(round((len(PUR_snps)/len(snp_id)*100), 2))+"% of all SNPs")

plot_snps(PUR_snps)

print()
print("Time at end:\t"+str(datetime.datetime.now().time()))
end = time.time()
print("Elapsed time: "+str(round((end-start)/60,2)))
print("----")
print()

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

chr22.close()
