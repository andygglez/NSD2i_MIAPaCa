
# geneBody plots
bash deepTools.plots.sh

# Fig 5E
python gazani_scripts/Fig5E.py --deseqDay5 deseq_tables/Day5.table.csv \
                --deseqDay9 deseq_tables/Day9.table.csv \
                --Veh_K27_D1_R1 bam_final/Vehicle_D1_Rep1_H3K27me3.clipped.bam \
                --Veh_K27_D1_R2 bam_final/Vehicle_D1_Rep2_H3K27me3.clipped.bam \
                --tpm_matrix tpm_symbol.txt --gtf references/hg38.gtf

# Fig 5F
python gazani_scripts/Fig5F.py --Veh_K27_D1_R1 bam_final/Vehicle_D1_Rep1_H3K27me3.clipped.bam \
                --Veh_K27_D1_R2 bam_final/Vehicle_D1_Rep2_H3K27me3.clipped.bam \
                --gtf references/hg38.gtf --deseqDay5 deseq_tables/Day5.table.csv \
                --deseqDay9 deseq_tables/Day9.table.csv --tpm_matrix tpm_symbol.txt

# Fig 5G
python gazani_scripts/Fig5G.py --enh_count references/genehancer_hg38_enhancer_count.txt \
                --deseqDay5 deseq_tables/Day5.table.csv --deseqDay9 deseq_tables/Day9.table.csv \
                --tpm_matrix tpm_symbol.txt
