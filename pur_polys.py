import os
import io
import ast
import yaml
import allel
import os.path
import itertools

import h5py
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import time        # to time the program
import datetime    # to print time at start of program

def plot_snps(PUR_snps, snp_id_cp, dest):

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
    plt.savefig(dest+'count_dist_chr.png')

    plt.figure()
    norm_dist = sns.distplot(pos[idx], bins=spatial_binnr, norm_hist=True)
    norm_dist.set_title("Normalized distribution")
    plt.savefig(dest+'norm_dist_chr.png')

    plt.figure()
    all_snps_count = sns.distplot(np.where(snp_id_cp[1,:] == 1), bins=spatial_binnr, kde=False, norm_hist=False)
    all_snps_count.set_title("Raw counts distribution of all SNPs")
    plt.savefig(dest+'raw_counts_snp_dist_chr.png')


print("Time at start:\t"+str(datetime.datetime.now().time()))

start = time.time()

with open('config.yaml', 'r') as f:
    content = yaml.load(f)

# file source and destination for i/o
SOURCE_PATH = content['source_dir']
DEST_PATH = content['output_dir']

# vcf source file name
VCF_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

# converted file name
HDF_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.h5'
# CSV_FILE = 'ALL.chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.csv'

# if the hdf5 file does not exist yet, start the conversion
if not os.path.exists(DEST_PATH+HDF_FILE):
    print("\nHDF file does not exist. Creating one...")
    allel.vcf_to_hdf5(SOURCE_PATH+VCF_FILE, DEST_PATH+HDF_FILE, fields='*', overwrite=True)

# load the data from the hdf5 file
print("\nLoading HDF file...")
chr = h5py.File(DEST_PATH+HDF_FILE, 'r')
chrom = 22

snp_id = chr['variants/ID'].value
z_arr = np.zeros((1, np.size(snp_id)))
snp_id_cp = np.vstack((snp_id, z_arr))     # add row full of zeros for future bool masking

pos = chr['variants/POS'].value
gt = chr['calldata/GT']
all_inds = chr['samples']

# read PUR individual list from tsv file (found on 1000Genomes)
# PUR_ids = np.genfromtxt(DEST_PATH+os.sep+"PUR_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
PUR_ids = np.genfromtxt(SOURCE_PATH+os.sep+"chr22_PUR_inds.txt", dtype=str, delimiter="\n")

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
        tmpPUR_snps.append(snp)
        snp_id_cp[1,i] = 1
    elif np.any(np.nansum(PUR_inds[i,:,:], axis = 0)) < 208 and np.any(np.nansum(nonPUR_inds[i,:,:], axis = 0)) == 2400:
        tmpPUR_snps.append(snp)
        snp_id_cp[1,i] = 1
    else:
        nonPUR_snps.append(snp)

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

with open(DEST_PATH+'chr'+str(chrom)+'PUR_SNPs.txt', 'w') as ff:
   ff.write(str(PUR_snps))

print()
print("PUR SNPs found:\t"+str(len(PUR_snps)))
print("Which is "+str(round((len(PUR_snps)/len(snp_id)*100), 2))+"% of all SNPs")

# plot_snps(PUR_snps)

print()
print("Time at end:\t"+str(datetime.datetime.now().time()))
end = time.time()
print("Elapsed time: "+str(round((end-start)/60,2)))
print("----")
print()

f.flush()
f.close()
ff.close()
chr.close()
