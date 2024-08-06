

marks=( H3K27Ac H3K27me3 H3K36me2 H3K36me3 H3K4me1 H3K4me3 H3K9me3 _NSD2 )
conditions=( NSD2i_D1 NSD2i_D5 NSD2i_D9 Vehicle_D1 Vehicle_D5 Vehicle_D9 )

bigWigs=$1

mkdir -p bigWigs_merged

parallel --jobs 32 --verbose bigwigCompare -b1 ${bigWigs}/{2}*Rep1*{1}* -b2 ${bigWigs}{2}*Rep2*{1}* --operation mean -p 6 -bs 10 -o bigWigs_merged/{2}_{1}.bigWig ::: ${marks[@]} ::: ${conditions[@]}