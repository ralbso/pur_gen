import yaml

with open('config.yaml', 'r') as f:
    content = yaml.load(f)

# file source and destination for i/o
SOURCE_PATH = content['source_dir']
DEST_PATH = content['output_dir']

HDF = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.h5'

def randomSample(dest_path, seed=301):
    import numpy as np
    import os
    from os import path
    import random

    pur_ids = np.genfromtxt(dest_path+"pur_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    ceu_ids = np.genfromtxt(dest_path+"ceu_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    gbr_ids = np.genfromtxt(dest_path+"gbr_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    ibs_ids = np.genfromtxt(dest_path+"ibs_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    pel_ids = np.genfromtxt(dest_path+"pel_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    tsi_ids = np.genfromtxt(dest_path+"tsi_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    yri_ids = np.genfromtxt(dest_path+"yri_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]

    random.seed(seed)
    ceu_sample = random.sample(list(ceu_ids),7)
    gbr_sample = random.sample(list(gbr_ids),7)
    ibs_sample = random.sample(list(ibs_ids),104)
    tsi_sample = random.sample(list(tsi_ids),12)
    pel_sample = random.sample(list(pel_ids),130)
    yri_sample = random.sample(list(yri_ids),130)

    if not os.path.exists(dest_path+'ceu_sample.txt'):
       with open(dest_path+'ceu_sample.txt', 'w') as ff:
           ff.write(str(np.array(ceu_sample)))

    if not os.path.exists(dest_path+'gbr_sample.txt'):
       with open(dest_path+'gbr_sample.txt', 'w') as ff:
           ff.write(str(np.array(gbr_sample)))

    if not os.path.exists(dest_path+'ibs_sample.txt'):
       with open(dest_path+'ibs_sample.txt', 'w') as ff:
           ff.write(str(ibs_sample))

    if not os.path.exists(dest_path+'pel_sampe.txt'):
       with open(dest_path+'pel_sample.txt', 'w') as ff:
           ff.write(str(pel_sample))

    if not os.path.exists(dest_path+'tsi_sample.txt'):
       with open(dest_path+'tsi_sample.txt', 'w') as ff:
           ff.write(str(tsi_sample))

    if not os.path.exists(dest_path+'yri_sample.txt'):
       with open(dest_path+'yri_sample.txt', 'w') as ff:
           ff.write(str(yri_sample))

    print(dest_path)

randomSample(DEST_PATH)

def createTextFiles(dest_path):
    import numpy as np
    import os
    from os import path

    populations = ["pur","ceu", "gbr","ibs","pel","tsi","yri"]
    pur_ids = np.genfromtxt(dest_path+"pur_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    ceu_ids = np.genfromtxt(dest_path+"ceu_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    gbr_ids = np.genfromtxt(dest_path+"gbr_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    ibs_ids = np.genfromtxt(dest_path+"ibs_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    pel_ids = np.genfromtxt(dest_path+"pel_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    tsi_ids = np.genfromtxt(dest_path+"tsi_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]
    yri_ids = np.genfromtxt(dest_path+"yri_igsr_samples.tsv", dtype=str, delimiter='\t')[1:,0]

    if not os.path.exists(dest_path+'pur_inds.txt'):
       with open(dest_path+'pur_inds.txt', 'w') as ff:
           ff.write(str(pur_ids))

    if not os.path.exists(dest_path+'ceu_inds.txt'):
       with open(dest_path+'ceu_inds.txt', 'w') as ff:
           ff.write(str(ceu_ids))

    if not os.path.exists(dest_path+'gbr_inds.txt'):
       with open(dest_path+'gbr_inds.txt', 'w') as ff:
           ff.write(str(gbr_ids))

    if not os.path.exists(dest_path+'ibs_inds.txt'):
       with open(dest_path+'ibs_inds.txt', 'w') as ff:
           ff.write(str(ibs_ids))

    if not os.path.exists(dest_path+'pel_inds.txt'):
       with open(dest_path+'pel_inds.txt', 'w') as ff:
           ff.write(str(pel_ids))

    if not os.path.exists(dest_path+'tsi_inds.txt'):
       with open(dest_path+'tsi_inds.txt', 'w') as ff:
           ff.write(str(tsi_ids))

    if not os.path.exists(dest_path+'yri_inds.txt'):
       with open(dest_path+'yri_inds.txt', 'w') as ff:
           ff.write(str(yri_ids))

createTextFiles(DEST_PATH)
