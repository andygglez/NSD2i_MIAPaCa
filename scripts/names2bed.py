import pandas as pd
import pyranges as pr
import argparse
import os

# # Parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('--geneNames', type=str, help='txt file with gene names (EnsembleID or gene symbols)')
parser.add_argument('--annotationGTF', type=str, help='Annotation file, GTF, to extract 1) the transcript coordinates and 2) the gene length quantiles')
parser.add_argument('--IDtype', type=str, help='Specify whether gene names are EnsemblID or gene symbols (ensembl | genesymbol)')
args = parser.parse_args()

gtf_folder="/".join(args.annotationGTF.split("/")[:-1]) if '/' in args.annotationGTF else "."
gtf_file=args.annotationGTF.split("/")[-1]

cache_file=gtf_folder+'/.cache.'+gtf_file+'.csv'
print("cache file should be:", cache_file)

if not os.path.exists(cache_file):
    
    # Reading GTF file
    gtf = pr.read_gtf(args.annotationGTF, as_df=True)

    # Filter GTF file by chromosomes
    chromosomes = ["chr"+str(i) for i in range(1,24)] ####### This will work for hg38 and mm10, not necesarilly with other genomes
    chromosomes.extend(['chrX', 'chrY'])

    # Take the transcripts from the GTF file
    transcriptsTable = gtf[['Chromosome','Feature','Start','End','Strand','gene_id','gene_name','gene_biotype','transcript_id','transcript_name','transcript_source','transcript_biotype','exon_id']]
    transcriptsTable = transcriptsTable.loc[transcriptsTable['Feature']=='transcript']
    transcriptsTable['Chromosome'] = [c for c in transcriptsTable['Chromosome']]
    transcriptsTable = transcriptsTable[['Chromosome', 'Start', 'End', 'gene_id', 'gene_name', 'transcript_id', 'Strand']]
    transcriptsTable.rename(columns={'Chromosome':'Chr', 'gene_id':'EnsemblID', 'gene_name':'GeneSymbol'}, inplace=True)
    transcriptsTable = transcriptsTable.loc[transcriptsTable['Chr'].isin(chromosomes)]
    transcriptsTable['GeneSymbol'] = transcriptsTable['GeneSymbol'].str.upper()
    transcriptsTable.reset_index(inplace=True, drop=True)

    transcriptsTable.to_csv(cache_file)
else:
    print("Reading existing cache file")
    transcriptsTable = pd.read_csv(cache_file)

# Define the column from the GTF file with IDs for the analysis (Ensembl IDs | GeneSymbols)
if args.IDtype == 'ensembl':
    column = 'EnsemblID'
elif args.IDtype == 'genesymbol': column = 'GeneSymbol'

# Read the txt file with gene names (Ensembl IDs | GeneSymbols)
with open(args.geneNames,'r') as file:
    genes = file.read()
    genes = genes.split('\n')
    genes = list(filter(lambda x: len(x)>1, genes))
genes = list(map(lambda x: x.upper(), genes))
input_total = len(genes)
######

# Take from the transcriptsTable only the genes in the txt file
df = transcriptsTable.loc[transcriptsTable[column].isin(genes)]
df.reset_index(inplace=True, drop=True)

notInGtf=list(set(genes).difference(set(transcriptsTable[column])))
nIG=pd.DataFrame()
nIG['abc']=notInGtf


# For each transcript df that belong to the same gene, select only one and change its Start and End values
# with the most upstream 'Start' and the most downstream 'End' of the transcripts of that gene
filtered = pd.DataFrame({c:[] for c in df.columns})
isoforms = list(df[column].unique())
for i in isoforms:
    df_ = df.loc[df[column]==i]
    df_.reset_index(inplace=True, drop=True)
    filtered.loc[filtered.shape[0]] = df_.loc[0]
    filtered.loc[filtered.shape[0]-1, 'Start'] = df_['Start'].min()
    filtered.loc[filtered.shape[0]-1, 'End'] = df_['End'].max()
    filtered.reset_index(inplace=True, drop=True)
filtered.reset_index(inplace=True, drop=True)

# Make a dataframe with the format of the bed file that will be used by deepTools
bed = filtered[['Chr','Start','End',column,'transcript_id','Strand']]
### Note: the 'transcript_id' column takes only the first 'transcript_id' that was present in the gtf file
### so it doesn't represent the real transcript taken for the experiment. The 'Start' and 'End' positions 
### correspond to the most upstream and downstream positions for all the transcripts of that gene

output_total = bed.shape[0]

name = args.geneNames.removesuffix('.txt')+'.bed'
bed.to_csv(name, header=False, index=False, sep='\t')
# nIG.to_csv('notInGTF.'+name, header=False, index=False, sep='\t')

percentage = output_total/input_total*100
print(f"Exported {percentage}% of genes from initial set.")