import sys
import os
import numpy as np
import pandas as pd

CC1=str(sys.argv[1])
CC2=str(sys.argv[2])
total_gene_file=str(sys.argv[3])
output_file = str(sys.argv[4])


total_gene = pd.read_csv(total_gene_file, sep='\t', header=0)
gff_file_CC1 = "/home/leeh7/iQTL/annotation_gft/gene_only/%s_anno_gene.txt" % CC1
gff_file_CC2 = "/home/leeh7/iQTL/annotation_gft/gene_only/%s_anno_gene.txt" % CC2
gff_CC1=pd.read_csv(gff_file_CC1, sep='\t', header=None, names=['chr','info','genes',CC1 + '_start',CC1+'_end','info1',CC1+'_strand',
                                                               'info2','gene_id'])
gff_CC2=pd.read_csv(gff_file_CC2, sep='\t', header=None, names=['chr','info','genes',CC2+'_start',CC2+'_end','info1',CC2+'_strand',
                                                               'info2','gene_id'])

gene_name_CC1=gff_CC1['gene_id'].str.partition(";Name=")
gene_name_CC2=gff_CC2['gene_id'].str.partition(";Name=")

gff_CC1['gene_name']=gene_name_CC1[2]
gff_CC2['gene_name']=gene_name_CC2[2]

gff_CC1['genes']= gene_name_CC1[0].str.partition("ID=gene:")[2]
gff_CC2['genes']= gene_name_CC2[0].str.partition("ID=gene:")[2]

CC1_gff=gff_CC1[['chr','genes',CC1 + '_start',CC1 + '_end',CC1 + '_strand','gene_name']]
CC2_gff=gff_CC2[['chr','genes',CC2 + '_start',CC2 + '_end',CC2 + '_strand','gene_name']]

total_gene = total_gene.rename(columns={'lg_multi': CC1 + '_multi', 'sm_multi': CC2 + '_multi'})

total_gene_counts=total_gene[['chr', CC1 + '_multi',CC2 + '_multi','mean_multi','genes',CC1 + '_unique', CC2 + '_unique','raw_total']]


merged_genes = CC1_gff.merge(CC2_gff, how='outer', on=['chr','genes','gene_name'])
merged_genes['chr']=merged_genes['chr'].astype(int)

merged_genes_counts = merged_genes.merge(total_gene_counts, how='inner', on=['chr','genes'])

column_to_move = merged_genes_counts.pop("gene_name")
merged_genes_counts.insert(2, "gene_name", column_to_move )

merged_genes_counts.to_csv(output_file, sep='\t', index=False, header=True)


