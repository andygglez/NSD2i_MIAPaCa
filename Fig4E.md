# convert aligned BAM files into 10kb binned BED files for each sample
```bash
bedtools bamtobed -i ${samp}.clipped.bam > reads.bed
sort -k1,1 -k2,2n reads.bed | gzip > ${samp}.bed.gz
bedtools intersect -a 10kb.bed -b ${samp}.bed.gz -sorted -c -nonamecheck > ${samp}.10kb.bed
```

# use R package ChIPbinner to normalize 10kb binned BED files and generate genic/intergenic scatterplot
```R
library(ChIPbinner)
# H3K36me2
## generate normalized bigWig files for each H3K36me2 replicate and sample (example below)
norm_bw(out_dir="outdir", genome_assembly = "hg38",chip_samp_binned_file = "MiaPaCa2.Vehicle_D1_Rep1.H3K36me2.10kb.bed",chip_samp_library_size = 20215104,use_control = FALSE,raw_count_cutoff = 10,pseudocount = 1)
norm_bw(out_dir="outdir", genome_assembly = "hg38",chip_samp_binned_file = "MiaPaCa2.Vehicle_D1_Rep2.H3K36me2.10kb.bed",chip_samp_library_size = 33407822,use_control = FALSE,raw_count_cutoff = 10,pseudocount = 1)
## merge replicates per condition (example below)
merge_norm_bw(are_R_objects = FALSE,
              rep1 = "MiaPaCa2.Vehicle_D1_Rep1.H3K36me2.10kb.norm.bw",
              rep2 = "MiaPaCa2.Vehicle_D1_Rep1.H3K36me2.10kb.norm.bw",
              merged_samples_label = "MiaPaCa2.Vehicle_D1.H3K36me2.10kb",
              out_dir = "outdir")
## genic/intergenic scatterplot
pow=0.3
h=4.25
w=3.85
max_x=100
max_y=100
# D1
genic_intergenic_scatterplot(
    out_dir = "outdir",
    genome_assembly = "hg38",
    cell_line = "MiaPaCa2",
    wildtype_samp_norm_bw = "MiaPaCa2.Vehicle_D1.H3K36me2.10kb.ms_norm.bw",
    treated_samp_norm_bw = "MiaPaCa2.NSD2i_D1.H3K36me2.10kb.ms_norm.bw",
    are_R_objects = FALSE,
    histone_mark = "H3K36me2",
    output_filename = "ms_norm.MiaPaCa2.Vehicle_D1.NSD2i_D1.H3K36me2.scatterplot",
    title_of_plot = "Day 1",
    xaxis_label = "DMSO",
    yaxis_label = "NSD2i",
    max_x = max_x,
    max_y = max_y,
    min_x = 0,
    min_y = 0,
    pow = pow,
    show_scales = TRUE,
    show_legend = FALSE,height_of_plot = h,width_of_plot = w
  )
# D5
genic_intergenic_scatterplot(
    out_dir = "outdir",
    genome_assembly = "hg38",
    cell_line = "MiaPaCa2",
    wildtype_samp_norm_bw = "MiaPaCa2.Vehicle_D5.H3K36me2.10kb.ms_norm.bw",
    treated_samp_norm_bw = "MiaPaCa2.NSD2i_D5.H3K36me2.10kb.ms_norm.bw",
    are_R_objects = FALSE,
    histone_mark = "H3K36me2",
    output_filename = "ms_norm.MiaPaCa2.Vehicle_D5.NSD2i_D5.H3K36me2.scatterplot",
    title_of_plot = "Day 5",
    xaxis_label = "DMSO",
    yaxis_label = "NSD2i",
    max_x = max_x,
    max_y = max_y,
    min_x = 0,
    min_y = 0,
    pow = pow,
    show_scales = TRUE,
    show_legend = FALSE,height_of_plot = h,width_of_plot = w
  )
# D9
genic_intergenic_scatterplot(
    out_dir = "outdir",
    genome_assembly = "hg38",
    cell_line = "MiaPaCa2",
    wildtype_samp_norm_bw = "MiaPaCa2.Vehicle_D9.H3K36me2.10kb.ms_norm.bw",
    treated_samp_norm_bw = "MiaPaCa2.NSD2i_D9.H3K36me2.10kb.ms_norm.bw",
    are_R_objects = FALSE,
    histone_mark = "H3K36me2",
    output_filename = "ms_norm.MiaPaCa2.Vehicle_D9.NSD2i_D9.H3K36me2.scatterplot",
    title_of_plot = "Day 9",
    xaxis_label = "DMSO",
    yaxis_label = "NSD2i",
    max_x = max_x,
    max_y = max_y,
    min_x = 0,
    min_y = 0,
    pow = pow,
    show_scales = TRUE,
    show_legend = FALSE,height_of_plot = h,width_of_plot = w
  )
# H3K27me3
## generate normalized bigWig files for each H3K36me2 replicate and sample (example below)
norm_bw(out_dir="outdir", genome_assembly = "hg38",chip_samp_binned_file = "MiaPaCa2.Vehicle_D1_Rep1.H3K27me3.10kb.bed",chip_samp_library_size = 31978009,use_control = FALSE,raw_count_cutoff = 10,pseudocount = 1)
norm_bw(out_dir="outdir", genome_assembly = "hg38",chip_samp_binned_file = "MiaPaCa2.Vehicle_D1_Rep2.H3K27me3.10kb.bed",chip_samp_library_size = 35411059,use_control = FALSE,raw_count_cutoff = 10,pseudocount = 1)
## merge replicates per condition (example below)
merge_norm_bw(are_R_objects = FALSE,
              rep1 = "MiaPaCa2.Vehicle_D1_Rep1.H3K27me3.10kb.norm.bw",
              rep2 = "MiaPaCa2.Vehicle_D1_Rep1.H3K27me3.10kb.norm.bw",
              merged_samples_label = "MiaPaCa2.Vehicle_D1.H3K27me3.10kb",
              out_dir = "outdir")
## genic/intergenic scatterplots
pow=0.36
# D1
genic_intergenic_scatterplot(
    out_dir = "outdir",
    genome_assembly = "hg38",
    cell_line = "MiaPaCa2",
    wildtype_samp_norm_bw = "MiaPaCa2.Vehicle_D1.H3K27me3.10kb.ms_norm.bw",
    treated_samp_norm_bw = "MiaPaCa2.NSD2i_D1.H3K27me3.10kb.ms_norm.bw",
    are_R_objects = FALSE,
    histone_mark = "H3K27me3",
    output_filename = "ms_norm.MiaPaCa2.Vehicle_D1.NSD2i_D1.H3K27me3.scatterplot",
    title_of_plot = "Day 1",
    xaxis_label = "DMSO",
    yaxis_label = "NSD2i",
    max_x = 50,
    max_y = 50,
    min_x = 0,
    min_y = 0,
    pow = pow,
    show_scales = TRUE,
    show_legend = FALSE,height_of_plot = 4,width_of_plot = 3.5
  )
# D5
genic_intergenic_scatterplot(
    out_dir = "outdir",
    genome_assembly = "hg38",
    cell_line = "MiaPaCa2",
    wildtype_samp_norm_bw = "MiaPaCa2.Vehicle_D5.H3K27me3.10kb.ms_norm.bw",
    treated_samp_norm_bw = "MiaPaCa2.NSD2i_D5.H3K27me3.10kb.ms_norm.bw",
    are_R_objects = FALSE,
    histone_mark = "H3K27me3",
    output_filename = "ms_norm.MiaPaCa2.Vehicle_D5.NSD2i_D5.H3K27me3.scatterplot",
    title_of_plot = "Day 5",
    xaxis_label = "DMSO",
    yaxis_label = "NSD2i",
    max_x = 50,
    max_y = 50,
    min_x = 0,
    min_y = 0,
    pow = pow,
    show_scales = TRUE,
    show_legend = FALSE,height_of_plot = 4,width_of_plot = 3.5
  )
# D9
genic_intergenic_scatterplot(
    out_dir = "outdir",
    genome_assembly = "hg38",
    cell_line = "MiaPaCa2",
    wildtype_samp_norm_bw = "MiaPaCa2.Vehicle_D9.H3K27me3.10kb.ms_norm.bw",
    treated_samp_norm_bw = "MiaPaCa2.NSD2i_D9.H3K27me3.10kb.ms_norm.bw",
    are_R_objects = FALSE,
    histone_mark = "H3K27me3",
    output_filename = "ms_norm.MiaPaCa2.Vehicle_D9.NSD2i_D9.H3K27me3.scatterplot",
    title_of_plot = "Day 9",
    xaxis_label = "DMSO",
    yaxis_label = "NSD2i",
    max_x = 50,
    max_y = 50,
    min_x = 0,
    min_y = 0,
    pow = pow,
    show_scales = TRUE,
    show_legend = FALSE,height_of_plot = 4,width_of_plot = 3.5
  )
```