import numpy as np
import pandas as pd

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

from hicrep import hicrepSCC
from hicrep.utils import readMcool

def get_scc_avg(fcool1, fcool2):
    
    if fcool1.endswith('mcool'):
        cool1, binSize1 = readMcool(fcool1, 100000)
    else:
        cool1, binSize1 = readMcool(fcool1, -1)
    if fcool2.endswith('mcool'):
        cool2, binSize2 = readMcool(fcool2, 100000)
    else:
        cool2, binSize2 = readMcool(fcool2, -1)

    for clr in cool1, cool2:
        first_20_chroms = clr.chroms()[:]['name'].values[:20]
        assert(np.all(first_20_chroms==np.array(['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
               'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               'chr16', 'chr17', 'chr18', 'chr19', 'chrX'])))  # confirm that the chromosomes are in order

    chromsizes = clr.chromsizes[:20].values
    chromsizes_norm = chromsizes/np.sum(chromsizes)  # sum to 1 (use as weights)
    
    # smoothing window half-size
    h = 1

    # maximal genomic distance to include in the calculation
    dBPMax = 500_000

    # whether to perform down-sampling or not 
    # if set True, it will bootstrap the data set # with larger contact counts to
    # the same number of contacts as in the other data set; otherwise, the contact 
    # matrices will be normalized by the respective total number of contacts
    bDownSample = False

    # compute the SCC score
    # this will result in a SCC score for each chromosome available in the data set
    # listed in the same order as the chromosomes are listed in the input Cooler files
    scc = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample)

    # only take numbered chromosomes and chrX
    scc = scc[:20]

    # average similarity score (weighted by chromosome size)
    scc_avg = np.sum(chromsizes_norm*scc)
    
    return scc_avg

### Get similarity matrix between 1B1/C36 Micro-C replicates

sample_names = ['1B1_rep1','1B1_rep2','1B1_rep6','1B1_rep7',
                'C36_rep1','C36_rep2','C36_rep6','C36_rep7',
                'G2_UT_rep1_genomewide','G2_UT_rep2_genomewide']
cooler_filenames = {}
for sample_name in sample_names:
    cooler_filenames[sample_name] = f'/mnt/coldstorage/jjusuf/MicroC_absloopquant_final/{sample_name}_all_100000.cool'

scc_matrix = pd.DataFrame(index=sample_names, columns=sample_names)
for i in range(len(sample_names)):
    for j in range(len(sample_names)):
        if j > i:
            print(f'working on {sample_names[i]} {sample_names[j]}: ', end='')
            scc_matrix.loc[sample_names[i], sample_names[j]] = get_scc_avg(cooler_filenames[sample_names[i]], cooler_filenames[sample_names[j]])
            print(f'{scc_matrix.loc[sample_names[i], sample_names[j]]}')

scc_matrix.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/MicroC_similarity_matrix_between_reps.csv')


### Get similarity matrix between Micro-C maps used for ultra-deep merged map

cooler_filenames_merged = {}
cooler_filenames_merged['1B1'] = '/mnt/coldstorage/jjusuf/MicroC_absloopquant_final/1B1_rep_merged_all_100000.cool'
cooler_filenames_merged['C36'] = '/mnt/coldstorage/jjusuf/MicroC_absloopquant_final/C36_rep_merged_all_100000.cool'
cooler_filenames_merged['C58'] = '/mnt/coldstorage/jjusuf/MicroC_absloopquant_final/C58_rep_merged_all_100000.cool'
cooler_filenames_merged['G2_UT_genomewide'] = '/mnt/coldstorage/jjusuf/MicroC_absloopquant_final/G2_UT_rep_merged_genomewide_all_100000.cool'
cooler_filenames_merged['Narducci2024'] = '/mnt/coldstorage/jjusuf/Past_MicroC_Experiments/all_WTgenome/mESC_Narducci2024_UNT_merged_updated_100000.cool'
cooler_filenames_merged['Paldi2024'] = '/mnt/coldstorage/jjusuf/Past_MicroC_Experiments/Paldi2024_DMSO/Paldi2024_DMSO_mC_rep_merged_100000.cool'
cooler_filenames_merged['Hsieh2022'] = '/mnt/coldstorage/jjusuf/Past_MicroC_Experiments/all_WTgenome/mESC_Hsieh2022_UT_merged.mcool'
cooler_filenames_merged['Hsieh2020'] = '/mnt/coldstorage/jjusuf/Past_MicroC_Experiments/all_WTgenome/mESC_Hsieh2020_WT_merged.mcool'
cooler_filenames_merged['Hansen2019'] = '/mnt/coldstorage/jjusuf/Past_MicroC_Experiments/all_WTgenome/mESC_Hansen2019_C59_WT_CTCF_all_merged.mcool'

sample_names_merged = list(cooler_filenames_merged.keys())

scc_matrix = pd.DataFrame(index=sample_names_merged, columns=sample_names_merged)

print('Working on similarity matrix for mega-merged')
for i in range(len(sample_names_merged)):
    for j in range(len(sample_names_merged)):
        if j > i:
            print(f'working on {sample_names_merged[i]} {sample_names_merged[j]}: ', end='')
            scc_matrix.loc[sample_names_merged[i], sample_names_merged[j]] = get_scc_avg(cooler_filenames_merged[sample_names_merged[i]], cooler_filenames_merged[sample_names_merged[j]])
            print(f'{scc_matrix.loc[sample_names_merged[i], sample_names_merged[j]]}')

scc_matrix.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/MicroC_similarity_matrix_mega_merged.csv')