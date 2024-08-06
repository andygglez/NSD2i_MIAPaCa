import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--tpm_matrix", type=str, help="tpm matrix")
parser.add_argument("--hg38_genes", type=str, help="bed file with coordinates of all protein coding genes")

args = parser.parse_args()

tpm = pd.read_table(args.tpm_matrix)
tpm['mean'] = tpm[['Day1_C1','Day1_C2','Day1_C3','Day1_C4']].mean(axis=1)
zero_names = list(tpm.loc[tpm['mean']==0]['gene'])

genes = pd.read_table(args.hg38_genes, names=['Chr', 'Start', 'End', 'GeneID', 'Transcript', 'Strand'])
genes = genes.loc[~genes['GeneID'].isin(zero_names)]
genes.to_csv("hg38.genes.filtExpression.bed", header=False, index=False, sep='\t')