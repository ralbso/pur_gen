def list_inds(chrom):
    import os
    import io
    import ast
    import os.path

    with open('D:/UPRM/GENOMES/chr'+str(chrom)+'_PUR_inds.txt', 'r') as ff:
       PUR_inds = [line.rstrip('\n') for line in ff]

    return(PUR_inds)
