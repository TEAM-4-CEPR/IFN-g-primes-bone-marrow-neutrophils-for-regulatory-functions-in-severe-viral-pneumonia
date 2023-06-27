#############################################
# Author : Ilango Guy                       #
# License : Team 4 CEPR                     #
# Corresponding Author : Christophe Paget   #
# Mail : christophe.paget@univ-tours.fr     #
#############################################

## Load Library
library(dplyr)
library("DESeq2")


## Load Dataset
df <- GSE109467_readCount_renamed_geneName

# Check name of column
names(df)

# Select only Preneu and Granulocyte-monocyte progenitors (GMP)
counts <- df %>% select(starts_with("PN") , starts_with("GMP"))
rownames(counts) <- df$X
counts


# Create Design Matrix
design <- data.frame(condition  =c("PN","PN","PN","GMP","GMP","GMP"))
rownames(design) <- colnames(counts %>% select(-gene_name))


#Checking matrix
all(rownames(design) %in% colnames(counts))
all(rownames(design) == colnames(counts))


## Run DEA
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = df3,
                              design = ~ condition)


featureData <- data.frame(gene=rownames(counts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

dds <- DESeq(dds)
res <- results(dds)
res <- results(dds, contrast=c("condition","PN","GMP"))
res

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_PN_vs_GMP", type="apeglm")
resLFC

# Filter DE Matrix
deg <- as.data.frame(resLFC) %>% filter(log2FoldChange != "NA")
deg <- deg %>% filter(padj <= 0.05)

# Add annotation
deg <- deg %>%
  mutate(annotation = case_when(
    (log2FoldChange > 0 & padj <= 0.05) ~ "Upregulated",
    (log2FoldChange < 0 & padj <= 0.05) ~ "Downregulated",
    (padj > 0.05) ~ "Non_sig",
  ))

deg$annotation <- as.factor(deg$annotation)
# VolcanoPlot
p<-ggplot(deg ,aes( x = log2FoldChange , y = -log10(padj) , color = annotation)) + geom_point()+theme_classic() + ggtitle("Volcano Plot GMP vs PN")+theme(plot.title = element_text(hjust = 0.5,face = "bold"))
p<- p+geom_vline(xintercept = 1.5, linetype= "dashed", color = "red", linewidth = 1)
p<-p+geom_vline(xintercept = -1.5, linetype= "dashed", color = "red", linewidth = 1)
p+ scale_color_manual(values = c("blue", "red"))

# Retrieve gene name within log2FC > 1.5 & log2FC < -1.5
gene1 <- rownames(deg %>% filter(log2FoldChange >= 1.5))
gene1
gene2 <- rownames(deg %>% filter(log2FoldChange <= -1.5))
gene <- c(gene1 , gene2)

check <-df %>% filter(X %in% gene)
write.csv(check , "gene_name.csv")




       