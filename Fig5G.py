
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.stats import mannwhitneyu

parser = argparse.ArgumentParser()
parser.add_argument("--enh_count", type=str, help="txt file with number of enhancers per gene")
parser.add_argument("--deseqDay5", type=str, help="Deseq2 report for day 5")
parser.add_argument("--deseqDay9", type=str, help="Deseq2 report for day 9")
parser.add_argument("--tpm_matrix", type=str, help="tpm matrix to extract controls")

args = parser.parse_args()


Day5 = pd.read_csv(args.deseqDay5)
Day5.rename(columns={'GeneSymbol':'Unnamed: 0', 'log2FoldChange_raw':'log2FoldChange', 'padj_raw':'padj'}, inplace=True)
Day5 = Day5.loc[Day5['padj']<0.05]

Day9 = pd.read_csv(args.deseqDay9)
Day9.rename(columns={'GeneSymbol':'Unnamed: 0', 'log2FoldChange_raw':'log2FoldChange', 'padj_raw':'padj'}, inplace=True)
Day9 = Day9.loc[Day9['padj']<0.05]


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
all_d1ctrl = all_d1ctrl.loc[all_d1ctrl.index.isin(Day5['Unnamed: 0'])]
all_d1ctrl.reset_index(inplace=True)

genes = pd.read_table(args.enh_count, names=['gene', 'count'])
genes = genes.loc[genes['gene'].isin(all_d1ctrl['gene'])]
genes.reset_index(inplace=True, drop=True)
enhDep = list(genes['gene'][:1000])

##########################
gt = pd.DataFrame()
gt['affected']=enhDep
gt.to_csv("Enhancer.Dependent.txt", header=False, index=False)
##########################

def find_closestControl(gene:str, tpm_matrixTest, tpm_matrixControl):
    gene_tpm = tpm_matrixTest.loc[tpm_matrixTest['gene']==gene, tpm_matrixTest.columns[1:]].values[0]
    distances = tpm_matrixControl.apply(lambda row: euclidean_distance(row[1:].values, gene_tpm), axis=1)
    return tpm_matrixControl.loc[distances.idxmin(), 'gene']


enhDep_tpm = all_d1ctrl.loc[all_d1ctrl['gene'].isin(enhDep)]
NotenhDep_tpm = all_d1ctrl.loc[all_d1ctrl['gene'].isin(genes['gene'][-2000:])]
controls = list(map(lambda x: find_closestControl(x, tpm_matrixTest=enhDep_tpm, tpm_matrixControl=NotenhDep_tpm), enhDep_tpm['gene']))

##########################
gt = pd.DataFrame()
gt['affected']=controls
gt.to_csv("Enhancer.Dependent.controls.txt", header=False, index=False)
##########################

controlLFC5 = Day5.loc[Day5['Unnamed: 0'].isin(controls)]['log2FoldChange']
enhDep_lfc5 = Day5.loc[Day5['Unnamed: 0'].isin(enhDep)]['log2FoldChange']

controlLFC9 = Day9.loc[Day9['Unnamed: 0'].isin(controls)]['log2FoldChange']
enhDep_lfc9 = Day9.loc[Day9['Unnamed: 0'].isin(enhDep)]['log2FoldChange']


# Creating plot
fig, ax = plt.subplots(figsize=(8, 4))

data=[controlLFC5, enhDep_lfc5, controlLFC9, enhDep_lfc9]

_, p_value5 = mannwhitneyu(controlLFC5, enhDep_lfc5)
_, p_value9 = mannwhitneyu(controlLFC9, enhDep_lfc9)

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

legend_labels=['Controls', 'Enhancer-dependent\ngenes']*2
legend_elements = [plt.Line2D([0], [0], color=color, lw=4, label=label) for color, label in zip(['tomato', 'blue'], legend_labels)]
ax.legend(handles=legend_elements, title='Gene Set', loc='center left', bbox_to_anchor=(0.8, 0.5))

# ax.set_xticks(np.arange(1,len(data)+1,1), x_labels)  # Set text labels.
ax.set_ylabel('LFC values in NSD2i')
ax.set_title(f"Effect of NSD2i on enhancer dependent genes")
ax.set_xlim(0.8,2.5)
ax.set_ylim(-4,3.1)
ax.set_xticks([1.15, 1.95], ['Day 5', 'Day9'])
ax.grid()
fig.savefig("images/Fig5G.svg")