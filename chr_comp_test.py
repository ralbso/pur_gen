import os
import allel
import os.path

import h5py
import numpy as np

PATH = 'C:'+os.sep+'Users'+os.sep+'rsoto'+os.sep+'Downloads'+os.sep
VCF_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
HDF_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.h5'
CSV_FILE = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.csv'

# VCF_FILE = 'ALL.chr22_GRCh38.genotypes.20170504.vcf.gz'
# HDF_FILE = 'ALL.chr22_GRCh38.genotypes.20170504.h5'
# CSV_FILE = 'ALL.chr22_GRCh38.genotypes.20170504.csv'

# allel.vcf_to_csv(PATH+VCF_FILE, PATH+CSV_FILE, fields='*')

if not os.path.exists(PATH+HDF_FILE):
    allel.vcf_to_hdf5(PATH+VCF_FILE, PATH+HDF_FILE, fields='*', overwrite=True)

chr22 = h5py.File(PATH+HDF_FILE, 'r')

snp_id = chr22['variants/ID']

afr_af = np.array(chr22['variants/AFR_AF'], dtype=np.float64)
amr_af = np.array(chr22['variants/AMR_AF'], dtype=np.float64)
sas_af = np.array(chr22['variants/SAS_AF'], dtype=np.float64)
eas_af = np.array(chr22['variants/EAS_AF'], dtype=np.float64)
eur_af = np.array(chr22['variants/EUR_AF'], dtype=np.float64)

sum_af = np.zeros([np.size(afr_af[:,1]),1], dtype=np.float64)
ind = np.zeros([np.size(afr_af[:,1]),1], dtype=np.float64)

# all_pops_af = np.concatenate((ind, afr_af, amr_af, sas_af, eas_af, eur_af), axis=1)
no_amr_af = np.concatenate((ind, afr_af, sas_af, eas_af, eur_af), axis=1)

thresh = 0

amr_snps = []
for i, af in enumerate(no_amr_af):
    if np.all(af == np.isnan):
        continue
    elif np.nansum(af) > thresh:
        no_amr_af[i,:] = np.nan
    elif np.nansum(af) <= thresh:
        amr_snps.append(i)

print(amr_snps)

rs_id = []
for i, abs_ind in enumerate(amr_snps):
    rs_id.append(snp_id[abs_ind])

rs_id.insert(0,len(rs_id)/np.size(afr_af[:,1]))

np.savetxt(str(PATH)+"AMR_SNPs_AF_LT_0.csv", rs_id, delimiter=",", fmt='%s', header='rs id')

chr22.close()
