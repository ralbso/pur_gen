import os
import io
import ast
import yaml
import allel
import os.path
import itertools

import h5py
import numpy as np
import numpy.ma as ma
import seaborn as sns
import matplotlib.pyplot as plt

with open('config.yaml', 'r') as f:
    content = yaml.load(f, Loader=yaml.FullLoader)

chrom = 22

# file source and destination for i/o
SOURCE_PATH = content['source_dir']
DEST_PATH = content['output_dir']

HDF_ALL_SNPS = 'PUR.chr22.phase3.20130502.genotypes.recode.h5'
HDF_PUR_SNPS = 'PUR.PUR_SNPs.chr22.phase3.20130502.genotypes.recode.h5'

# if the hdf5 file does not exist yet, start the conversion
if not os.path.exists(DEST_PATH+HDF_ALL_SNPS):
    print("\nHDF file does not exist. Creating one...")
    VCF_ALL_SNPS = 'PUR.chr22.phase3.20130502.genotypes.recode.vcf'
    allel.vcf_to_hdf5(DEST_PATH+VCF_ALL_SNPS, DEST_PATH+HDF_ALL_SNPS, fields='*', overwrite=True)

# if the hdf5 file does not exist yet, start the conversion
if not os.path.exists(DEST_PATH+HDF_PUR_SNPS):
    print("\nHDF file does not exist. Creating one...")
    VCF_PUR_SNPS = 'PUR.PUR_SNPs.chr22.phase3.20130502.genotypes.recode.vcf'
    allel.vcf_to_hdf5(DEST_PATH+VCF_PUR_SNPS, DEST_PATH+HDF_PUR_SNPS, fields='*', overwrite=True)

all_snps = h5py.File(DEST_PATH+HDF_ALL_SNPS, 'r')
pur_snps = h5py.File(DEST_PATH+HDF_PUR_SNPS, 'r')

snp_id_all = all_snps['variants/ID'].value
pos_all = all_snps['variants/POS'].value

snp_id_pur = pur_snps['variants/ID'].value
pos_pur = pur_snps['variants/POS'].value

snp_mask = np.in1d(snp_id_all, snp_id_pur, assume_unique=True)
print("ALL")
print(snp_id_all[:10])
print()
print("PUR")
print(snp_id_pur[:10])
print()
print("Mask")
print(snp_mask[:10])
