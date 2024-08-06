import pybedtools
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--deseqDay5", type=str, help="Deseq2 report for day 5")
parser.add_argument("--deseqDay9", type=str, help="Deseq2 report for day 9")
parser.add_argument("--Veh_K27_D1_R1", type=str, help="Bam file with sample from K27me3 on day 1, rep1")
parser.add_argument("--Veh_K27_D1_R2", type=str, help="Bam file with sample from K27me3 on day 1, rep2")
parser.add_argument("--gtf", type=str, help="GTF annotation file to extract gene coordinates")
parser.add_argument("--tpm_matrix", type=str, help="tpm matrix to extract controls")

args = parser.parse_args()

Day5 = pd.read_csv(args.deseqDay5)
Day5.rename(columns={'GeneSymbol':'Unnamed: 0', 'log2FoldChange_raw':'log2FoldChange', 'padj_raw':'padj'}, inplace=True)

Day9 = pd.read_csv(args.deseqDay9)
Day9.rename(columns={'GeneSymbol':'Unnamed: 0', 'log2FoldChange_raw':'log2FoldChange', 'padj_raw':'padj'}, inplace=True)

if 'references' not in os.listdir():
    os.mkdir('references')
if 'images' not in os.listdir():
    os.mkdir('images')

print("Calling peaks for H3K27me3 Day 1")
## I'm calling peaks on the parental sample of K27me3, so I can intersect the peaks with promoters and identify the PRC2 target genes 
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

print("---------------------------------------------")

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

promoters = pd.concat([pos_strand, neg_strand])

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

subprocess.run("rm fc.matrix*", shell=True, check=True)

print("Taking the first 500 promoters with more H3K27me3")
tpm_matrixK27me3 = tpm_matrixK27me3.loc[tpm_matrixK27me3['Geneid'].isin(set(Day5.loc[Day5['padj']<0.05]['Unnamed: 0']).intersection(set(Day9.loc[Day9['padj']<0.05]['Unnamed: 0'])))]

prc2_list = list(tpm_matrixK27me3['Geneid'][:500]) # Select first 1000 genes
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
all_d1ctrl = all_d1ctrl.loc[all_d1ctrl.index.isin(set(Day5.loc[Day5['padj']<0.05]['Unnamed: 0']).intersection(set(Day9.loc[Day9['padj']<0.05]['Unnamed: 0'])))]
all_d1ctrl.reset_index(inplace=True)

hg38_genes = pybedtools.BedTool.from_dataframe(hg38_genes)

def find_closestControl(gene:str, tpm_matrixTest, tpm_matrixControl):
    gene_tpm = tpm_matrixTest.loc[tpm_matrixTest['gene']==gene, tpm_matrixTest.columns[1:]].values[0]
    distances = tpm_matrixControl.apply(lambda row: euclidean_distance(row[1:].values, gene_tpm), axis=1)
    return tpm_matrixControl.loc[distances.idxmin(), 'gene']

def get_DataX_PRC2(deseq_result, tpm_matrix, input_genes_lfc=None, mode=1):
    from scipy.stats import mannwhitneyu

    deseq_filt = deseq_result.loc[deseq_result['padj']<0.05]
    
    prc2_lfc = deseq_filt.loc[deseq_filt['Unnamed: 0'].isin(prc2_list)]['log2FoldChange'] ##################################################

    if mode==1:
        tpm_matrix_filt = tpm_matrix.loc[tpm_matrix['gene'].isin(deseq_filt['Unnamed: 0'])]        # tpm_matrix_filt = tpm_matrix_filt.loc[tpm_matrix_filt[tpm_matrix_filt.columns[1:]].mean(axis=1)>0]

        prc2Targets_tpmTable = tpm_matrix_filt.loc[tpm_matrix_filt['gene'].isin(prc2_list)] ###################333
        Notprc2Targets_tpmTable = tpm_matrix_filt.loc[~ tpm_matrix_filt['gene'].isin(prc2_list)] ################333
        
        controls = list(map(lambda x: find_closestControl(x, tpm_matrixTest=prc2Targets_tpmTable, tpm_matrixControl=Notprc2Targets_tpmTable), prc2Targets_tpmTable['gene']))
        controlLFC = deseq_filt.loc[deseq_filt['Unnamed: 0'].isin(controls)]['log2FoldChange']
        
        return controlLFC, prc2_lfc, controls

    if mode==2:

        deseq_filt = deseq_result.loc[deseq_result['padj']<0.05]

        prc2_lfc = deseq_filt.loc[deseq_filt['Unnamed: 0'].isin(prc2_list)]['log2FoldChange'] 
        controlLFC = deseq_filt.loc[deseq_filt['Unnamed: 0'].isin(input_genes_lfc)]['log2FoldChange']

        print(f"Received: {len(input_genes_lfc)}; Output: {len(controlLFC)}")
        return controlLFC, prc2_lfc, input_genes_lfc


controlLFC_PRC2_Day5, PRC2_Day5, control_names = get_DataX_PRC2(deseq_result=Day5, tpm_matrix=all_d1ctrl, mode=1)

####################################################################################
# y_labels=['PRC2 targets', 'Controls']
# data = [PRC2_Day5, controlLFC_PRC2_Day5]
# _, pvalue = mannwhitneyu(PRC2_Day5, controlLFC_PRC2_Day5)

# # Creating plot
# fig, ax = plt.subplots(figsize=(8, 4))

# print('-'*50)
# print(f"{len(data[0])} : {len(data[1])}")


# # Create a list of colors for the boxplots based on the number of features you have
# boxplots_colors = ['blue', 'tomato','forestgreen', 'plum'] if len(data)==4 else ['blue', 'tomato']

# # Boxplot data
# box_width=0.2
# bp = ax.boxplot(data, patch_artist = True, vert = False, widths=[box_width for _ in range(len(data))])

# # Change to the desired color and add transparency
# for patch, color in zip(bp['boxes'], boxplots_colors):
#     patch.set_facecolor(color)
#     patch.set_alpha(0.6)

# # Create a list of colors for the violin plots based on the number of features you have
# violin_colors = ['blue', 'tomato','forestgreen', 'plum'] if len(data)==4 else ['blue', 'tomato']

# # Violinplot data
# vp = ax.violinplot(data, points=500, 
#             showmeans=False, showextrema=False, showmedians=False, vert=False)

# for idx, b in enumerate(vp['bodies']):
#     # Get the center of the plot
#     m = np.mean(b.get_paths()[0].vertices[:, 0])
#     # Modify it so we only see the upper half of the violin plot
#     b.get_paths()[0].vertices[:, 1] = np.clip(b.get_paths()[0].vertices[:, 1], idx+1, idx+2)
#     # Change to the desired color
#     b.set_color(violin_colors[idx])

# # Create a list of colors for the scatter plots based on the number of features you have
# scatter_colors = ['blue', 'tomato','forestgreen', 'plum'] if len(data)==4 else ['blue', 'tomato']

# # Scatterplot data
# for idx, features in enumerate(data):
#     # Add jitter effect so the features do not overlap on the y-axis
#     y = np.full(len(features), idx + .8)
#     idxs = np.arange(len(y))
#     out = y.astype(float)
#     out.flat[idxs] += np.random.uniform(low=-.05, high=.05, size=len(idxs))
#     y = out
#     ax.scatter(features, y, s=.3, c=scatter_colors[idx])

# ax.set_yticks(np.arange(1,len(data)+1,1), y_labels)  # Set text labels.
# ax.set_xlabel('LFC values in NSD2i')
# ax.set_title(f"Effect of NSD2i on PRC2 targets")
# if pvalue:
#     ax.text(-1, 1.5 ,f'p-value: {pvalue:.3e} (Mann-Whitney U test)')

# bbox_style = {'facecolor': 'snow', 'alpha': 0.5, 'boxstyle': "round, pad=0.5", 'ec': 'black'}
# ax.text(2, 1, f'mean: {np.mean(data[0]):.2f} \nmedian: {np.median(data[0]):.2f}', bbox = bbox_style)
# ax.text(2, 2, f'mean: {np.mean(data[1]):.2f} \nmedian: {np.median(data[1]):.2f}', bbox = bbox_style)

# ax.set_xlim(-2.5,3.2)
# ax.grid()
# fig.savefig("images/Fig5E.svg", bbox_inches='tight')


####################################################################################

controlLFC_PRC2_Day9, PRC2_Day9, names_control_lfc = get_DataX_PRC2(deseq_result=Day9, input_genes_lfc=control_names, tpm_matrix=all_d1ctrl, mode=2)

# Creating plot
fig, ax = plt.subplots(figsize=(8, 4))

data=[controlLFC_PRC2_Day5, PRC2_Day5, controlLFC_PRC2_Day9, PRC2_Day9]

print(f"controlLFC_PRC2_Day5 {len(controlLFC_PRC2_Day5)}")
print(f"PRC2_Day5 {len(PRC2_Day5)}")

print(f"controlLFC_PRC2_Day9 {len(controlLFC_PRC2_Day9)}")
print(f"PRC2_Day9 {len(PRC2_Day9)}")

_, p_value5 = mannwhitneyu(controlLFC_PRC2_Day5, PRC2_Day5)
_, p_value9 = mannwhitneyu(controlLFC_PRC2_Day9, PRC2_Day9)

print('-'*50)
if len(data)==2:
    print(f"{len(data[0])} : {len(data[1])}")
if len(data)==4:
    print(f"{len(data[0])} : {len(data[2])} ------- {len(data[1])} : {len(data[3])}")

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

legend_labels=['Controls', 'PRC2 targets']*2
legend_elements = [plt.Line2D([0], [0], color=color, lw=4, label=label) for color, label in zip(['tomato', 'blue'], legend_labels)]
ax.legend(handles=legend_elements, title='Gene Set', bbox_to_anchor=(1, 0.1))

# ax.set_xticks(np.arange(1,len(data)+1,1), x_labels)  # Set text labels.
ax.set_ylabel('LFC values in NSD2i')
ax.set_title(f"Effect of NSD2i on PRC2 targets")
ax.set_xlim(0.8,2.5)
ax.set_ylim(-4,3.1)
ax.set_xticks([1.15, 1.95], ['Day 5', 'Day9'])
ax.grid()
fig.savefig("images/Fig5E.svg", bbox_inches='tight')