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

with open('config.yaml', 'r') as f:
    content = yaml.load(f)

# file source and destination for i/o
SOURCE_PATH = content['source_dir']
DEST_PATH = content['output_dir']

VCF_FILE = 'PUR.PUR_SNPs.chr22.phase3.20130502.genotypes.recode.vcf'
HDF_FILE = 'PUR.PUR_SNPs.chr22.phase3.20130502.genotypes.recode.h5'

# if the hdf5 file does not exist yet, start the conversion
if not os.path.exists(DEST_PATH+HDF_FILE):
    print("\nHDF file does not exist. Creating one...")
    allel.vcf_to_hdf5(DEST_PATH+VCF_FILE, DEST_PATH+HDF_FILE, fields='*', overwrite=True)

data = h5py.File(DEST_PATH+HDF_FILE, 'r')

snp_id = data['variants/ID'].value
pos = data['variants/POS'].value
