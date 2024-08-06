marks=( H3K27me3 H3K36me2 H3K36me3 H3K4me3 H3K4me1 H3K9me3 H3K27Ac _NSD2 )

mkdir -p deepTools

declare -A ylims=( [_NSD2]="0 0.2"
                   [H3K27me3]="0 11"
                   [H3K36me2]="0 50"
                   [H3K36me3]="0 4.2"
                   [H3K4me1]="0 15"
                   [H3K4me3]="0 1.5"
                   [H3K9me3]="0 6"
                   [H3K27Ac]="0 0.2"
)

declare -A heatmaplims=( [_NSD2]="0 0.32"
                   [H3K27me3]="0 38"
                   [H3K36me2]="0 110"
                   [H3K36me3]="0 18"
                   [H3K4me1]="0 55"
                   [H3K4me3]="0 3.5"
                   [H3K9me3]="0 40"
                   [H3K27Ac]="0 0.6"
)


plot () {
    local m=$1

    bwVehicle=$(ls MS_normalized/Vehicle*${m}*bigWig)
    bwNSD2i=$(ls MS_normalized/NSD2*${m}*bigWig)


    all_bw=( ${bwVehicle[@]} ${bwNSD2i[@]} )

    
    IFS=" " read -r -a lims <<< "${ylims[${m}]}"
    
    IFS=" " read -r -a heatlims <<< "${heatmaplims[${m}]}"
    
    

############################### Vehicle

############# Without labels references/hg38.protein.coding.genes.bed
    computeMatrix scale-regions -R hg38.genes.filtExpression.bed \
        -S ${bwVehicle[@]} \
        --binSize 250 -a 4500 -b 4500 -m 6000 \
        --missingDataAsZero --skipZeros \
        -bl ~/references/hg38/hg38.blacklist.bed \
        -p "max/2" --samplesLabel "" "" "" \
        -o deepTools/matrix.veh.${m}.gz
    

    plotHeatmap -m deepTools/matrix.veh.${m}.gz \
        -o deepTools/heatmap.veh.${m}.svg \
        --dpi 600 \
        --colorMap 'Reds' --perGroup --colorList 'white,red,blue' \
        -T "Vehicle" --yMin ${lims[0]} --yMax ${lims[1]} \
        --zMin ${heatlims[0]} --zMax ${heatlims[1]} \
        --regionsLabel ${m}
        
############# With labels
        
    computeMatrix scale-regions -R hg38.genes.filtExpression.bed \
        -S ${bwVehicle[@]} \
        --binSize 250 -a 0 -b 0 -m 6000 \
        --missingDataAsZero --skipZeros \
        -bl ~/references/hg38/hg38.blacklist.bed \
        -p "max/2" --samplesLabel "D1" "D5" "D9" \
        -o deepTools/matrix.veh.centered0.${m}.gz

        
############################### NSD2i

############# Without labels

    computeMatrix scale-regions -R hg38.genes.filtExpression.bed \
        -S ${bwNSD2i[@]} \
        --binSize 250 -a 4500 -b 4500 -m 6000 \
        --missingDataAsZero --skipZeros \
        -bl ~/references/hg38/hg38.blacklist.bed \
        -p "max/2" --samplesLabel "" "" "" \
        -o deepTools/matrix.nsd2i.${m}.gz

    plotHeatmap -m deepTools/matrix.nsd2i.${m}.gz \
        -o deepTools/heatmap.nsd2i.${m}.svg \
        --dpi 600 \
        --colorMap 'Reds' --perGroup --colorList 'white,red,blue' \
        -T "NSD2i"  --yMin ${lims[0]} --yMax ${lims[1]} \
        --zMin ${heatlims[0]} --zMax ${heatlims[1]} \
        --regionsLabel ${m}

############# With labels

    computeMatrix scale-regions -R hg38.genes.filtExpression.bed \
        -S ${bwNSD2i[@]} \
        --binSize 250 -a 0 -b 0 -m 6000 \
        --missingDataAsZero --skipZeros \
        -bl ~/references/hg38/hg38.blacklist.bed \
        -p "max/2" --samplesLabel "D1" "D5" "D9" \
        -o deepTools/matrix.nsd2i.centered0.${m}.gz


}

for m in ${marks[@]}; do

    plot ${m}

done