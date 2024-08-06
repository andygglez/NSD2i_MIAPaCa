marks=( H3K27me3 H3K36me2 )

mkdir -p deepTools

declare -A ylims=([H3K27me3]="0 12" [H3K36me2]="0 25")

declare -A heatmaplims=( [H3K27me3]="0 36" [H3K36me2]="0 75")

plot () {

    local m=$1

    bwVehicle=$(ls MS_normalized/Vehicle*${m}*bigWig)
    bwNSD2i=$(ls MS_normalized/NSD2*${m}*bigWig)


    all_bw=( ${bwVehicle[@]} ${bwNSD2i[@]} )

    echo ${all_bw[@]}
    echo "-----------------------------------"


############################### Vehicle

    computeMatrix scale-regions -R hg38.igr_20kb.bed \
        -S ${bwVehicle[@]} \
        --binSize 1000 -m 20000 -a 20000 -b 20000 \
        --missingDataAsZero --skipZeros \
        -bl ~/references/hg38/hg38.blacklist.bed \
        -p "max/2" --samplesLabel "" "" "" \
        -o deepTools/matrix.veh.intergenic.${m}.gz
    
    IFS=" " read -r -a lims <<< ${ylims[${m}]}
    
    IFS=" " read -r -a heatlims <<< "${heatmaplims[${m}]}"
    

    plotHeatmap -m deepTools/matrix.veh.intergenic.${m}.gz \
        -o deepTools/heatmap.veh.intergenic.${m}.svg \
        --dpi 600 \
        --colorMap 'Reds' --perGroup --colorList 'white,red,blue' --startLabel "" --endLabel "" \
        -T "Vehicle, intergenic regions" --yMin ${lims[0]} --yMax ${lims[1]} \
        --zMin ${heatlims[0]} --zMax ${heatlims[1]} \
        --regionsLabel ${m}
        
############################### NSD2
    
     computeMatrix scale-regions -R hg38.igr_20kb.bed \
        -S ${bwNSD2i[@]} \
        --binSize 1000 -m 20000 -a 20000 -b 20000 \
        --missingDataAsZero --skipZeros \
        -bl ~/references/hg38/hg38.blacklist.bed \
        -p "max/2" --samplesLabel "" "" "" \
        -o deepTools/matrix.nsd2i.intergenic.${m}.gz

    plotHeatmap -m deepTools/matrix.nsd2i.intergenic.${m}.gz \
        -o deepTools/heatmap.nsd2i.intergenic.${m}.svg \
        --dpi 600 \
        --colorMap 'Reds' --perGroup --colorList 'white,red,blue' --startLabel "" --endLabel "" \
        -T "NSD2i, intergenic regions" --yMin ${lims[0]} --yMax ${lims[1]} \
        --zMin ${heatlims[0]} --zMax ${heatlims[1]} \
        --regionsLabel ${m}

}

for m in ${marks[@]}; do

    plot ${m}

done
