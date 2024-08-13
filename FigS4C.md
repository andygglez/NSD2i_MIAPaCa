# run featureCounts for IDR-filtered NiceSeq peaks
```bash
# create SAF files for pooled peaks
cat idr/*.bed | sort -k1,1 -k2,2n | bedtools merge | tee MiaPaCa2.merged_idr_peaks.bed | wc -l
#transform merged.bed to SAF format for featureCounts
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' MiaPaCa2.merged_idr_peaks.bed > MiaPaCa2.merged_idr_peaks.saf
# run featureCounts
featureCounts -a MiaPaCa2.merged_idr_peaks.saf -F SAF -o MiaPaCa2.NiceSeq.counts.txt -T 6 -p MiaPaCa2_Vehicle_Day_1_NicE-Seq__S1_L001.trimmed.srt.nodup.bam MiaPaCa2_Vehicle_Day_1_NicE-Seq__S2_L001.trimmed.srt.nodup.bam MiaPaCa2_NSD2i_Day_1_NicE-Seq__S3_L001.trimmed.srt.nodup.bam /MiaPaCa2_NSD2i_Day_1_NicE-Seq__S4_L001.trimmed.srt.nodup.bam MiaPaCa2_Vehicle_Day_5_NicE-Seq__S5_L001.trimmed.srt.nodup.bam MiaPaCa2_Vehicle_Day_5_NicE-Seq__S6_L001.trimmed.srt.nodup.bam MiaPaCa2_NSD2i_Day_5_NicE-Seq__S7_L001.trimmed.srt.nodup.bam MiaPaCa2_NSD2i_Day_5_NicE-Seq__S8_L001.trimmed.srt.nodup.bam MiaPaCa2_Vehicle_Day_9_NicE-Seq__S9_L001.trimmed.srt.nodup.bam MiaPaCa2_Vehicle_Day_9_NicE-Seq__S10_L001.trimmed.srt.nodup.bam MiaPaCa2_NSD2i_Day_9_NicE-Seq__S11_L001.trimmed.srt.nodup.bam MiaPaCa2_NSD2i_Day_9_NicE-Seq__S12_L001.trimmed.srt.nodup.bam
```

# run DESeq2 on raw counts and generate PCA
```R
df <- fread("MiaPaCa2.NiceSeq.counts.txt") %>% as.data.frame() %>% column_to_rownames("Geneid")
columns_to_drop <- c('Chr', 'Start', 'End','Strand','Length')
df <- df[, -which(names(df) %in% columns_to_drop)]
# rename samples
colnames(df) <- gsub("_L001.trimmed.srt.nodup.bam|23CB2[0-9][0-9]_|_NicE-Seq_|MiaPaCa2_","",colnames(df))
colnames(df)
# drop outlier sample
df <- df[!grepl('NSD2i_Day_9_S12', names(df))]
df <- df[!grepl('Day_9', names(df))]
# create matrix
mat <- df %>% as.matrix()
mdat <- data.frame(kind=colnames(mat)) %>% `rownames<-` (colnames(mat))
mdat$condition <- mdat$kind
mdat$condition <- gsub("_S1[0-9]|_S[0-9]","",mdat$condition)
# run DESEq2
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = mdat,
                              design = ~ condition)
dds <- DESeq(dds)
# generate PCA data
keep <- rowSums(counts(dds)) >= 10
dds.filt <- dds[keep,]
vsd <- vst(dds.filt)
pca_dt <- DESeq2::plotPCA(vsd, intgroup = "condition", ntop = 500,returnData = T)
DESeq2::plotPCA(vsd, intgroup = "condition", ntop = 500)
# add shapes to differentiate group
pca_dt$group <- "NSD2i"
pca_dt[1:2,3] <- "DMSO"
pca_dt[5:6,3] <- "DMSO"
pca_dt$timepoint <- "Day 1"
pca_dt[5:8,6] <- "Day 5"
# plot PCA
pca_dt %>%
  ggplot(aes(x = PC1, y = PC2, color = timepoint,label=name,shape=group)) +
  geom_point(size = 3, show.legend = T) +
  geom_hline(yintercept = 0,alpha=0.1) +
  geom_vline(xintercept = 0,alpha=0.1) +
  xlab(label = "PC1 (40%)") +
  ylab(label = "PC2 (23%)") +
  labs(title = "Top 500 most variable Nice-Seq peaks",colour="") +
  scale_shape_manual(values=c(16, 17),name="Treatment",labels=c("DMSO","NSD2i")) +
  scale_color_manual(values = c("blue","forestgreen","red"),name="Timepoint") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=12,family="Helvetica",colour = "black"),
    axis.text.y=element_text(size=12,family="Helvetica",colour = "black"),
    axis.text.x=element_text(size=12,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "black",size=12,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
        size = 12, color = "black",family = "Helvetica"),
    legend.text=element_text(size=12,family = "Helvetica",color = "black"),
    legend.title=element_text(size=12,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    legend.key.size = unit(0.1, 'cm'),
    plot.margin = unit(c(1.2,1.2,1.2,1.2), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0.1,'cm'),
    panel.spacing.x = unit(0.1,'cm'))
ggsave(filename = "NiceSeq_PCA_top500.pdf",device = "pdf",units = "cm",width = 15,height=12,dpi = 600,bg="white")
```