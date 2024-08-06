# NSD2i_MIAPaCa
A repository for scripts from Majewski's lab in our contribution to NSD2i paper

In each of our scripts we provide commented code to produce our figures in the paper. The folder **references/** contains files that we used in some of our analyses, such as:
- MS_gozani.yaml: MS coefficients used to normalize the bigWig files
- hg38.igr_20kb.bed: bed file with intergenic regions used for creating deepTools plots
- hg38.protein.coding.genes.bed: bed file with genic regions used for creating deepTools plots

The folder **scripts/** provides some accesory scripts that are called from the main ones producing the images.

The **run.sh** file specifies how to call each of the scripts to produce the images
