import pybedtools
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--Veh_K27_D1_R1", type=str, help="Bam file with sample from K27me3 on day 1, rep1")
parser.add_argument("--Veh_K27_D1_R2", type=str, help="Bam file with sample from K27me3 on day 1, rep2")
parser.add_argument("--gtf", type=str, help="GTF annotation file to extract gene coordinates")
parser.add_argument("--deseqDay5", type=str, help="Deseq2 report for day 5")
parser.add_argument("--deseqDay9", type=str, help="Deseq2 report for day 9")
parser.add_argument("--tpm_matrix", type=str, help="tpm matrix to extract controls")

args = parser.parse_args()


if 'references' not in os.listdir():
    os.mkdir('references')
if 'images' not in os.listdir():
    os.mkdir('images')

print("Calling peaks for H3K27me3 Day 1")
# I'm calling peaks on the parental sample of K27me3, so I can intersect the peaks with promoters and identify the PRC2 target genes 
call_peaks1 = subprocess.run(f"macs2 callpeak -t {args.Veh_K27_D1_R1} -f BAM -g hs --keep-dup all --outdir . --name K27me3.rep1", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
call_peaks2 = subprocess.run(f"macs2 callpeak -t {args.Veh_K27_D1_R2} -f BAM -g hs --keep-dup all --outdir . --name K27me3.rep2", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

awk = subprocess.run("awk 'BEGIN {OFS=\"\t\"} {print $1,$2,$3}' K27me3.rep1_peaks.narrowPeak > K27me3.rep1_peaks.bed", shell=True, check=True)
awk = subprocess.run("awk 'BEGIN {OFS=\"\t\"} {print $1,$2,$3}' K27me3.rep2_peaks.narrowPeak > K27me3.rep2_peaks.bed", shell=True, check=True)

cat = subprocess.run("cat K27me3.rep1_peaks.bed >> K27me3.merged.bed", shell=True, check=True)
cat = subprocess.run("cat K27me3.rep2_peaks.bed >> K27me3.merged.bed", shell=True, check=True)

bed_sort = subprocess.run("bedtools sort -i K27me3.merged.bed > K27me3.merged.sorted.bed", shell=True, check=True)
bed_merge = subprocess.run("bedtools merge -i K27me3.merged.sorted.bed > K27me3.merged.sorted.final.bed", shell=True, check=True)

mv = subprocess.run("mv K27me3.merged.sorted.final.bed references/", shell=True, check=True)
rm = subprocess.run("rm K27me3*rep* K27me3.merged.bed K27me3.merged.sorted.bed", shell=True, check=True)

call_peaks1_broad = subprocess.run(f"macs2 callpeak --broad -t {args.Veh_K27_D1_R1} -f BAM -g hs --keep-dup all --outdir . --name K27me3.broad.rep1", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
call_peaks2_broad = subprocess.run(f"macs2 callpeak --broad -t {args.Veh_K27_D1_R2} -f BAM -g hs --keep-dup all --outdir . --name K27me3.broad.rep2", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

awk = subprocess.run("awk 'BEGIN {OFS=\"\t\"} {print $1,$2,$3}' K27me3.broad.rep1_peaks.broadPeak > K27me3.broad.rep1_peaks.bed", shell=True, check=True)
awk = subprocess.run("awk 'BEGIN {OFS=\"\t\"} {print $1,$2,$3}' K27me3.broad.rep2_peaks.broadPeak > K27me3.broad.rep2_peaks.bed", shell=True, check=True)

cat = subprocess.run("cat K27me3.broad.rep1_peaks.bed >> K27me3.broad.merged.bed", shell=True, check=True)
cat = subprocess.run("cat K27me3.broad.rep2_peaks.bed >> K27me3.broad.merged.bed", shell=True, check=True)

bed_sort = subprocess.run("bedtools sort -i K27me3.broad.merged.bed > K27me3.broad.merged.sorted.bed", shell=True, check=True)
bed_merge = subprocess.run("bedtools merge -i K27me3.broad.merged.sorted.bed > K27me3.broad.merged.sorted.final.bed", shell=True, check=True)

mv = subprocess.run("mv K27me3.broad.merged.sorted.final.bed references/", shell=True, check=True)
rm = subprocess.run("rm K27me3*rep* K27me3.broad.merged.bed K27me3.broad.merged.sorted.bed", shell=True, check=True)

cat = subprocess.run("cat references/K27me3.merged.sorted.final.bed > references/K27me3.broad.narrow.merged.bed", shell=True, check=True)
cat = subprocess.run("cat references/K27me3.broad.merged.sorted.final.bed >> references/K27me3.broad.narrow.merged.bed", shell=True, check=True)

bed_sort = subprocess.run("bedtools sort -i references/K27me3.broad.narrow.merged.bed | bedtools merge -d 30000 > references/K27me3.broadNarrow.merged.distance.bed", shell=True, check=True)
rm = subprocess.run("rm references/K27me3.broad.narrow.merged.bed references/K27me3.broad.merged.sorted.final.bed references/K27me3.merged.sorted.final.bed", shell=True, check=True)

subprocess.run("awk \'$3==\"gene\" {print $0}\' "+args.gtf+" | awk \'$14==\"\\\"protein_coding\\\";\" {print $0}\' | awk \'{print $18}\' > references/hg38.protein.coding.genes.txt", shell=True, check=True)
subprocess.run("sed -i \'s/;//g\' references/hg38.protein.coding.genes.txt", shell=True, check=True)
subprocess.run("sed -i \'s/\"//g\' references/hg38.protein.coding.genes.txt", shell=True, check=True)

subprocess.run("python scripts/names2bed.py --geneNames references/hg38.protein.coding.genes.txt --annotationGTF "+args.gtf+" --IDtype genesymbol", shell=True, check=True)
subprocess.run("rm references/hg38.protein.coding.genes.txt", shell=True, check=True)


hg38_genes = pd.read_table("references/hg38.protein.coding.genes.bed", names=['Chr', 'Start', 'End', 'GeneID', 'TranscriptID', 'Strand'])
hg38_genes = hg38_genes[['Chr','Start','End','Strand','GeneID']]
# pos_strand = table.loc[table['Strand']=='+'][['Chr', 'Start']]
pos_strand = hg38_genes.loc[hg38_genes['Strand']=='+'][['Chr', 'Start', 'End', 'Strand', 'GeneID']]
pos_strand['PromStart']=pos_strand['Start']-5000
pos_strand['PromEnd']=pos_strand['Start']+800

# # neg_strand = table.loc[table['Strand']=='-'][['Chr', 'End']]
neg_strand = hg38_genes.loc[hg38_genes['Strand']=='-'][['Chr', 'Start', 'End', 'Strand', 'GeneID']]
neg_strand['PromStart']=neg_strand['End']-800
neg_strand['PromEnd']=neg_strand['End']+5000



hg38_genes = pd.read_table("references/hg38.protein.coding.genes.bed", names=['Chr', 'Start', 'End', 'GeneID', 'TranscriptID', 'Strand'])
hg38_genes = hg38_genes[['Chr','Start','End','Strand','GeneID']]
# pos_strand = table.loc[table['Strand']=='+'][['Chr', 'Start']]
pos_strand = hg38_genes.loc[hg38_genes['Strand']=='+'][['Chr', 'Start', 'End', 'Strand', 'GeneID']]
pos_strand['PromStart']=pos_strand['Start']-5000
pos_strand['PromEnd']=pos_strand['Start']+800

# # neg_strand = table.loc[table['Strand']=='-'][['Chr', 'End']]
neg_strand = hg38_genes.loc[hg38_genes['Strand']=='-'][['Chr', 'Start', 'End', 'Strand', 'GeneID']]
neg_strand['PromStart']=neg_strand['End']-800
neg_strand['PromEnd']=neg_strand['End']+5000


promoters = pd.concat([pos_strand, neg_strand])

Day5 = pd.read_csv(args.deseqDay5)
Day5.rename(columns={'GeneSymbol':'Unnamed: 0', 'log2FoldChange_raw':'log2FoldChange', 'padj_raw':'padj'}, inplace=True)

Day9 = pd.read_csv(args.deseqDay9)
Day9.rename(columns={'GeneSymbol':'Unnamed: 0', 'log2FoldChange_raw':'log2FoldChange', 'padj_raw':'padj'}, inplace=True)
k27me3 = pybedtools.BedTool("references/K27me3.broadNarrow.merged.distance.bed")

promoters = pybedtools.BedTool.from_dataframe(promoters)
PRC2_targets_promoters = promoters.intersect(k27me3, wa=True)
PRC2_targets_promoters_DF = PRC2_targets_promoters.to_dataframe()
PRC2_targets_promoters_DF.rename(columns={'name':'strand', 'score':'name'}, inplace=True)


subprocess.run("awk \'BEGIN {OFS=\"\t\"} {print $5,$1,$2,$3, \".\"}\' references/hg38.prc2.targets.promoters.bed > hg38.prc2.targets.promoters.saf", shell=True, check=True)
subprocess.run(f"featureCounts -a hg38.prc2.targets.promoters.saf -F SAF -o fc.matrix -T 2 -p {args.Veh_K27_D1_R1} {args.Veh_K27_D1_R2}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
subprocess.run("rm hg38.prc2.targets.promoters.saf", shell=True, check=True)
subprocess.run("awk \'NR>1\' fc.matrix > fc.matrix.filt", shell=True, check=True)


fC_matrix = pd.read_table("fc.matrix.filt")
len_dict = dict(zip(fC_matrix['Geneid'], fC_matrix['Length']))



def tpm(fc_matrix, len_dict):
    col1 = fc_matrix.columns[0]
    # Normalize by gene length
    fc = fc_matrix.copy()
    for i in range(len(fc)):
        if len_dict.get(fc.loc[i][col1]):
            fc.loc[i, fc.columns[1:]] = fc.loc[i, fc.columns[1:]]/len_dict[fc.loc[i][col1]]
        else: fc.loc[i, fc.columns[1:]] = np.nan
    sums = [fc[c].sum()/1000000 for c in fc.columns[1:]]
    for c in range(1, len(fc.columns)):
        fc[fc.columns[c]] = fc[fc.columns[c]]/sums[c-1]
    return fc

tpm_matrixK27me3 = tpm(fC_matrix[['Geneid', args.Veh_K27_D1_R1, args.Veh_K27_D1_R2]], len_dict=len_dict)
tpm_matrixK27me3['mean'] = tpm_matrixK27me3[[args.Veh_K27_D1_R1, args.Veh_K27_D1_R2]].mean(axis=1)
tpm_matrixK27me3.sort_values(by='mean', ascending=False, inplace=True)



tpm_matrixK27me3 = tpm_matrixK27me3.loc[tpm_matrixK27me3['Geneid'].isin(set(Day5.loc[Day5['padj']<0.05]['Unnamed: 0']).intersection(set(Day9.loc[Day9['padj']<0.05]['Unnamed: 0'])))]

print(f"Taking the first 500 promoters with more H3K27me3 {tpm_matrixK27me3.shape[0]}")
prc2_list = list(tpm_matrixK27me3['Geneid'][:500]) # Select first 500 genes
prc2_targets_coord = hg38_genes.loc[hg38_genes['GeneID'].isin(prc2_list)]
PRC2_targets = list(prc2_targets_coord['GeneID'])
PRC2_coordinates = hg38_genes.loc[hg38_genes['GeneID'].isin(PRC2_targets)]
PRC2_coordinates.rename(columns={'GeneID':'ID'}, inplace=True)


# Read and modify the tpm matrices
print("Reading and modifying tpm matrices")
tpm = pd.read_table(args.tpm_matrix, header=0, index_col=0)
tpm = tpm.groupby(tpm.index).mean()

## log2tpm
tpm_normed = np.log2(np.divide(*(tpm.iloc[:,np.where(tpm.columns.str.contains('_N'))[0]]+1).\
                    align((tpm.iloc[:,np.where(tpm.columns.str.contains('_C'))[0]]+1), axis=0)))

tpm_normed.columns = tpm_normed.columns.str.replace('N', 'Rep').str.replace('Day', 'D') + '_log2tpm'

## merged gene set
gene_uni = set(tpm_normed.index)

## slecting random geens from control d1
def euclidean_distance(row1, row2):
    return np.sqrt(np.sum((row1 - row2) ** 2))

# base expression
all_d1ctrl = tpm.loc[gene_uni,tpm.columns.str.contains('Day1_C')]
all_d1ctrl.reset_index(inplace=True)

hg38_genes = pybedtools.BedTool.from_dataframe(hg38_genes)


def expand_bed(df_bed:pd.DataFrame, dist:int=10000)->pd.DataFrame:
    # This function WILL include the gene itself
    chromosomes = ["chr"+str(i) for i in range(1,24)]
    chromosomes.extend(['chrX','chrY'])

    # Extension
    new = df_bed.copy()
    new['Start2'] = [int(i)-dist if int(i)-dist>0 else 0 for i in new['Start']]
    new['End2'] = [int(i)+dist for i in new['End']]  
    new = new.loc[new['Chr'].isin(chromosomes)]

    return new[['Chr', 'Start2', 'End2', 'ID']]

def find_closestControl(gene:str, tpm_matrixTest, tpm_matrixControl):
    gene_tpm = tpm_matrixTest.loc[tpm_matrixTest['gene']==gene, tpm_matrixTest.columns[1:]].values[0]
    distances = tpm_matrixControl.apply(lambda row: euclidean_distance(row[1:].values, gene_tpm), axis=1)
    return tpm_matrixControl.loc[distances.idxmin(), 'gene']

def get_DataX(dist, deseq_result, tpm_matrix, input_genes_lfc=None, mode=1):
    from scipy.stats import mannwhitneyu

    # Take the genes and expand the flanks
    genes_expanded = expand_bed(PRC2_coordinates, dist=dist)
    genes_expanded = pybedtools.BedTool.from_dataframe(genes_expanded)

    # Take the genes that are touched by the expansion, including PRC2 targets
    genes_expanded2 = hg38_genes.intersect(genes_expanded, wa=True) # This is going to include PRC2 target genes
    PRC2_genes = pybedtools.BedTool.from_dataframe(PRC2_coordinates)
    genes_touched = genes_expanded2.intersect(PRC2_genes, v=True) # Take the genes affected by the expansion of K27me3, exclude PRC2 targets

    genes_touched = genes_touched.to_dataframe().drop_duplicates()
    genes_touched.rename(columns={'name':'strand','score':'name'}, inplace=True)

    deseq_filt = deseq_result.loc[deseq_result['padj']<0.05]

    #######################################################################################################################################################
    
    tpm_matrix = tpm_matrix.loc[tpm_matrix['gene'].isin(deseq_filt['Unnamed: 0'])]

    # Take only the genes within the filtered tpm matrix and deseq_filtered
    genes_touched = genes_touched.loc[(genes_touched['name'].isin(tpm_matrix['gene']))]


    # Take the genes located in the flanks of PRC2 targets and find their logFoldChange values
    genesTouched_list = genes_touched['name']
        
    allhg38_genes = hg38_genes.to_dataframe()
    allhg38_genes.rename(columns={'name':'strand','score':'name'}, inplace=True)

    # Genes from the controls cannot be in the test set nor in the list of PRC2 targets
    genesNotTouched_list = list(allhg38_genes.loc[(~ allhg38_genes['name'].isin(genesTouched_list)) & (~ allhg38_genes['name'].isin(PRC2_targets)) & (allhg38_genes['name'].isin(tpm_matrix['gene']))]['name'])


    genesTouched_lfc = deseq_filt.loc[deseq_filt['Unnamed: 0'].isin(genesTouched_list)]['log2FoldChange']

    # Finding the controls
    # tpm_filtered = tpm_matrix.loc[(tpm_matrix['id'].isin(genesTouched_list)) | (tpm_matrix['id'].isin(genesNotTouched_list))] # Take only the tpm of genes that have lfc to report

    genesTouched_tpmTable = tpm_matrix.loc[tpm_matrix['gene'].isin(genesTouched_list)]
    genesNotTouched_tpmTable = tpm_matrix.loc[tpm_matrix['gene'].isin(genesNotTouched_list)]
    
    controls = list(map(lambda x: find_closestControl(x, tpm_matrixTest=genesTouched_tpmTable, tpm_matrixControl=genesNotTouched_tpmTable), genesTouched_tpmTable['gene']))
    
    ##########################
    gt = pd.DataFrame()
    gt['affected']=genesTouched_tpmTable['gene']
    gt.to_csv("PRC2.neighbors.txt", header=False, index=False)
    ##########################

    ##########################
    gt = pd.DataFrame()
    gt['affected']=controls
    gt.to_csv("PRC2.neighbors.controls.txt", header=False, index=False)
    ##########################

    controlLFC = deseq_filt.loc[deseq_filt['Unnamed: 0'].isin(controls)]['log2FoldChange']
    
    _, p_value = mannwhitneyu(controlLFC, genesTouched_lfc)

    return controlLFC, genesTouched_lfc, p_value

print("Obtaining plots")
distance=400000
controls, affected, pvalue = get_DataX(dist=distance, deseq_result=Day5, tpm_matrix=all_d1ctrl)
controls9, affected9, pvalue9 = get_DataX(dist=distance, deseq_result=Day9, tpm_matrix=all_d1ctrl)

print(f"pv-5: {pvalue}, pv-9: {pvalue9}")
print(len(prc2_list))

# Creating plot
fig, ax = plt.subplots(figsize=(8, 4))

data=[controls, affected, controls9, affected9]

_, p_value5 = mannwhitneyu(controls, affected)
_, p_value9 = mannwhitneyu(controls9, affected9)

# Create a list of colors for the boxplots based on the number of features you have
boxplots_colors = ['tomato', 'blue']*2

# Boxplot data
box_width=0.1
positions=[1,1.3, 1.8, 2.1]
# bp = ax.boxplot(data, patch_artist = True, vert = False, widths=[box_width for _ in range(len(data))])
bp = ax.boxplot(data, patch_artist = True, vert = True, widths=[box_width for _ in range(len(data))], positions=positions)

# Change to the desired color and add transparency
for patch, color in zip(bp['boxes'], boxplots_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

# Create a list of colors for the violin plots based on the number of features you have
violin_colors = ['tomato', 'blue']*2

# Violinplot data
# vp = ax.violinplot(data, points=500, 
#             showmeans=False, showextrema=False, showmedians=False, vert=False)
vp = ax.violinplot(data, points=500, 
            showmeans=False, showextrema=False, showmedians=False, vert=True, positions=positions, widths=[box_width*2.5 for _ in range(len(data))])

for idx, b in enumerate(vp['bodies']):
    # Get the center of the plot
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # Modify it so we only see the upper half of the violin plot
    # b.get_paths()[0].vertices[:, 1] = np.clip(b.get_paths()[0].vertices[:, 1], idx+1, idx+2)
    # Change to the desired color
    b.set_color(violin_colors[idx])

# Create a list of colors for the scatter plots based on the number of features you have
scatter_colors = ['tomato', 'blue']*2

bbox_style = {'facecolor': 'snow', 'alpha': 0.5, 'boxstyle': "round, pad=0.5", 'ec': 'black'}

ax.text(np.mean(positions)-0.52, 2, f'p-value: {p_value5:.2e}', bbox = bbox_style)
ax.text(np.mean(positions)+0.6, 2, f'p-value: {p_value9:.2e}', bbox = bbox_style)

for i in range(len(data)):
    ax.text(positions[i]+0.05, -3, f'mean: {np.mean(data[i]):.2f} \nmedian: {np.median(data[i]):.2f}', bbox = bbox_style)

legend_labels=['Controls', 'PRC2 targets\nneighbors']*2
legend_elements = [plt.Line2D([0], [0], color=color, lw=4, label=label) for color, label in zip(['tomato', 'blue'], legend_labels)]
ax.legend(handles=legend_elements, title='Gene Set', loc='center left', bbox_to_anchor=(0.8, 0.5))

# ax.set_xticks(np.arange(1,len(data)+1,1), x_labels)  # Set text labels.
ax.set_ylabel('LFC values in NSD2i')
ax.set_title(f"Effect of NSD2i on genes neighboring PRC2 targets. Distance: 400 kb")
ax.set_xlim(0.8,2.5)
ax.set_ylim(-4,3.1)
ax.set_xticks([1.15, 1.95], ['Day 5', 'Day9'])
ax.grid()
fig.savefig("images/Fig5F.svg", bbox_inches='tight')