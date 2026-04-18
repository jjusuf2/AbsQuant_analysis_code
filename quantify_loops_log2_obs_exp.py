import numpy as np
import pandas as pd
import cooler
from tqdm import tqdm
from multiprocess import Pool

## setup ##

# load the list of filtered loops
filtered_loops_df = pd.read_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/filtered_loops.csv', index_col=0)

# load the coolers & P(s) curves

sample_names = ['1B1_rep1','1B1_rep2','1B1_rep6','1B1_rep7',
                'C36_rep1','C36_rep2','C36_rep6','C36_rep7']
coolers = {}
for sample_name in sample_names:
    coolers[sample_name] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/MicroC_absloopquant_final/{sample_name}_all_1000.cool')

def log2oe_score(clr, chrom, left, right):
    try:
        window = clr.matrix().fetch(f'{chrom}:{left-20000}-{left+20000}',f'{chrom}:{right-20000}-{right+20000}')
    except:
        return np.nan
    exp1 = window[:5, :5]
    exp2 = window[-5:,-5:]
    exp_signal = (np.mean(exp1)+np.mean(exp2))/2
    obs = window[18:23,18:23]
    obs_signal = np.mean(obs)
    return np.log2(obs_signal / exp_signal)

log2oe_scores = pd.DataFrame(index=filtered_loops_df.index, columns=sample_names)

for sample_name in sample_names:
    print(f'working on {sample_name}')
    clr = coolers[sample_name]
    
    for i in tqdm(range(len(filtered_loops_df))):
        chrom, left, right, size = filtered_loops_df.iloc[i]
        score = log2oe_score(clr, chrom, left, right)
        log2oe_scores.loc[filtered_loops_df.index[i], sample_name] = score

        if i % 1000 == 0:
            log2oe_scores.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/log2_obs_over_exp_scores.csv')

    log2oe_scores.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/log2_obs_over_exp_scores.csv')
log2oe_scores.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/log2_obs_over_exp_scores.csv')