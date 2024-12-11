from statistics import mean
import numpy as np
import pandas as pd
import sys
import os

lg_bed=str(sys.argv[2])
sm_bed=str(sys.argv[3])
uniq_bed=str(sys.argv[1])
CC1=str(sys.argv[4])
CC2=str(sys.argv[5])
mouse=str(sys.argv[6])
output_path=str(sys.argv[7])

unique_bed_file = pd.read_csv(uniq_bed, sep='\t', header=None, names=['chr','start','stop','sense','genes','counts'])
lg_multi = pd.read_csv(lg_bed, sep='\t', header=None, names=['chr','start','stop','sense','genes','lg_multi'])
sm_multi = pd.read_csv(sm_bed, sep='\t', header=None, names=['chr','start','stop','sense','genes','sm_multi'])

lg_multi = lg_multi.groupby(['chr','genes']).lg_multi.sum().reset_index()
sm_multi = sm_multi.groupby(['chr','genes']).sm_multi.sum().reset_index()
sm_multi['genes'] = sm_multi['genes'].str.replace(CC2, '').astype(str)
lg_multi['genes'] = lg_multi['genes'].str.replace(CC1, '').astype(str)
lg_multi['chr'] = lg_multi['chr'].str.replace(r'chr', '').astype(str)
sm_multi['chr'] = sm_multi['chr'].str.replace(r'chr', '').astype(str)

full_current_lib = lg_multi.merge(sm_multi, how='inner', on=['chr','genes'])
full_current_lib['mean_multi'] = full_current_lib[['lg_multi','sm_multi']].mean(axis=1)

unique_bed_file = unique_bed_file.groupby(['chr','genes']).counts.sum().reset_index()
# grab LG-specific unique counts
lg_unique = unique_bed_file[unique_bed_file['chr'].str.contains(CC1 + 'chr')].rename(columns={'counts' : CC1 + '_unique'})
lg_unique['chr'] = lg_unique['chr'].str.replace(CC1 + r'chr', '') ##.astype(int)
lg_unique['genes'] = lg_unique['genes'].str.replace(CC1, '').astype(str)
# grab SM-specific counts
sm_unique = unique_bed_file[unique_bed_file['chr'].str.contains(CC2 + 'chr')].rename(columns={'counts' : CC2 + '_unique'})
sm_unique['chr'] = sm_unique['chr'].str.replace(CC2 + r'chr', '') ##.astype(int)
sm_unique['genes'] = sm_unique['genes'].str.replace(CC2, '').astype(str)

full_current_lib = full_current_lib.merge(lg_unique, how='inner', on=['chr','genes'])
full_current_lib = full_current_lib.merge(sm_unique, how='inner', on=['chr','genes'])

full_current_lib['raw_total'] = full_current_lib['mean_multi'] + full_current_lib[CC1 + '_unique'] + full_current_lib[CC2 + '_unique'] 
# filter out genes that are not expressed (0 total read counts)
full_current_lib = full_current_lib[full_current_lib['raw_total'] > 0]

output_filename = output_path + mouse + "_total_counts.txt" 
full_current_lib.to_csv(output_filename, sep='\t', index=False, header=True)
