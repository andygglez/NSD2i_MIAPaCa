# Specify folder with bigWig files. In this case: bigWigs/
bash scripts/merge.bw.sh bigWigs/

mkdir -p MS_normalized/
mkdir -p bigWigs_merged/

parallel --jobs 32 --verbose Rscript scripts/ms_norm_bw.R --bigWig {1} --MS_yaml references/MS_gozani.yaml --removeSuffix .bigWig --exportFolder MS_normalized ::: $(ls bigWigs_merged/*{H3K27me3,H3K36me2,H3K36me3,H3K4me1,H3K9me3}*)

# Copy non-MS-normalized bigWigs to the final folder
cp bigWigs_merged/*{_NSD2,K4me3,K27Ac}*bigWig MS_normalized/

# bash scripts/deepTools.genic.sh
bash scripts/deepTools.intergenic.sh
