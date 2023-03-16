library(ggplot2)
library(dplyr)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(data.table)
library(ggpubr)
library(viridis)
library(scales)
library(ggbeeswarm)
library(gmodels)
library(DESeq2)

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
library(rlang)
library(clusterProfiler)
library(enrichplot)
font_import()

#organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

library(MetBrewer)
library(msigdbr)
library(ggh4x)

library(stringr)

palette <- met.brewer(name="Veronese", n=7, type="discrete")

met.brewer(name="Veronese", n=7, type="discrete")

palette2 <- c(palette[3], palette[7])


theme_science <- function (base_size = 12, base_family = "Arial Black") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}

makeStars <- function(x){
  stars <- c("***", "**", "*", "ns")
  vec <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

scaleFUN <- function(x) sprintf("%.2f", x)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbviridis <- c("#440154FF", "#31688EFF", "#35B779FF", "#E69F00")

setwd("/Users/cottr/Box Sync/APEX-ADAR/RNAseq/")


#read data
df_r <- fread("all.gene_counts.tsv")

length(unique(df_r$ensembl_gene_id))
length(unique(df_r$entrezgene))
length(unique(df_r$external_gene_name))

df_r$id <- paste(row.names(df_r), df_r$ensembl_gene_id, sep = "_")

df_r_idents <- df_r[,c("id", "external_gene_name", "ensembl_gene_id")]

#read ATF4 targets from _______________________________
atf4_targets <- read.delim("elife-42940-supp1-v1.txt")

atf4_targets <- data.frame(V1 = toupper(atf4_targets[1:69,2]))

atf4_targets <- merge(atf4_targets, df_r_idents, by.x = "V1", by.y = "external_gene_name")

atf4_targets <- atf4_targets$ensembl_gene_id

#read core_isgs 
Core_ISGs <- read.delim("Core_ISGs.txt")

Core_ISGs <- merge(Core_ISGs, df_r_idents, by.x = "ISGs", by.y = "external_gene_name")

ISGS <- as.character(Core_ISGs$ensembl_gene_id)

h_isgs_a <- read.delim("HALLMARK_INTERFERON_ALPHA_RESPONSE.v2022.1.Hs.grp", header = FALSE)
h_isgs_g <- read.delim("HALLMARK_INTERFERON_GAMMA_RESPONSE.v2022.1.Hs.grp", header = FALSE)

h_isgs <- unique(rbind(h_isgs_a, h_isgs_g))

nfkb <- read.delim("HALLMARK_TNFA_SIGNALING_VIA_NFKB.v2022.1.Hs.grp", header = FALSE)
nfkb <- subset(nfkb, !V1 %in% h_isgs$V1)

nfkb <- merge(nfkb, df_r_idents, by.x = "V1", by.y = "external_gene_name")

nfkb <- as.character(nfkb$ensembl_gene_id)

upr <- read.delim("HALLMARK_UNFOLDED_PROTEIN_RESPONSE.v2022.1.Hs.grp", header = FALSE)
upr <- merge(upr, df_r_idents, by.x = "V1", by.y = "external_gene_name")

upr <- as.character(upr$ensembl_gene_id)




######################################################
######
###### INTERACTION in SKBR3
######
######
######################################################


data_skbr3 <- as.matrix(df_r[,8:15])

row.names(data_skbr3) <- df_r$id

head(data_skbr3)

coldata_skbr3 <- DataFrame(shrna1 = factor(c("a", "a", "s", "s","a", "a", "s", "s")),
                          shrna2 = factor(c("d","s","d","s","d","s","d","s")),
                          row.names=as.character(colnames(data_skbr3)))
coldata_skbr3

ddsHTSeq_skbr3 <- DESeqDataSetFromMatrix(countData= data_skbr3, 
                                        colData = coldata_skbr3, 
                                        design = ~ shrna1 + shrna2 + shrna1:shrna2)
ddsHTSeq_skbr3

ddsHTSeq_skbr3$shrna1
ddsHTSeq_skbr3$shrna1 = relevel(ddsHTSeq_skbr3$shrna1, "s")
ddsHTSeq_skbr3$shrna1

ddsHTSeq_skbr3$shrna2
ddsHTSeq_skbr3$shrna2 = relevel(ddsHTSeq_skbr3$shrna2, "s")
ddsHTSeq_skbr3$shrna2

dds_skbr3 <- DESeq(ddsHTSeq_skbr3)

resultsNames(dds_skbr3)

keep <- rowSums(counts(dds_skbr3)) >= 10
dds_skbr3 <- dds_skbr3[keep,]


interaction_skbr3 <- results(dds_skbr3, name = "shrna1a.shrna2d")

interaction_skbr3_shrink <- lfcShrink(dds_skbr3, type="apeglm", coef = 4)

plotMA(interaction_skbr3)
plotMA(interaction_skbr3_shrink)

sum(interaction_skbr3$padj <0.01, na.rm=TRUE)
sum(interaction_skbr3_shrink$padj <0.01, na.rm=TRUE)

interaction_skbr3_df <- as.data.frame(interaction_skbr3)
interaction_skbr3_shrink_df <- as.data.frame(interaction_skbr3_shrink)

interaction_skbr3_df$id<- row.names(interaction_skbr3_df)
interaction_skbr3_shrink_df$id<- row.names(interaction_skbr3_shrink_df)

interaction_skbr3_df <- merge(interaction_skbr3_df, df_r, by = "id")
interaction_skbr3_shrink_df <- merge(interaction_skbr3_shrink_df, df_r, by = "id")

write.csv(interaction_skbr3_shrink_df, "interaction_skbr3_shrink_df.csv")


interaction_skbr3_shrink_df$Sig <- ifelse(interaction_skbr3_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(interaction_skbr3_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_interaction.tiff", height=4, width=4, unit="in")




interaction_skbr3_shrink_df$ISG <- grepl(paste(ISGS, collapse = "|"), interaction_skbr3_shrink_df$ensembl_gene_id)
interaction_skbr3_shrink_df$ATF4_targets <- grepl(paste(atf4_targets, collapse = "|"), interaction_skbr3_shrink_df$ensembl_gene_id)

ggplot(interaction_skbr3_shrink_df, aes(log2FoldChange, -log10(padj), colour = ISG, alpha = ISG)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_interaction_ISG.tiff", height=4, width=4, unit="in")


ggplot(interaction_skbr3_shrink_df, aes(log2FoldChange, -log10(padj), colour = ATF4_targets, alpha = ATF4_targets)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_interaction_ATF4.tiff", height=4, width=4, unit="in")


interaction_skbr3_shrink_df$rank <- -log10(interaction_skbr3_shrink_df$padj)*sign(interaction_skbr3_shrink_df$log2FoldChange)

interaction_skbr3_shrink_df <- do.call(data.frame, lapply(interaction_skbr3_shrink_df, function(x) replace(x, is.infinite(x), NA)))


# we want the log2 fold change 
gene_list_interaction_skbr3 <- interaction_skbr3_shrink_df$rank

# name the vector
names(gene_list_interaction_skbr3) <- interaction_skbr3_shrink_df$ensembl_gene_id

# omit any NA values 
gene_list_interaction_skbr3<-na.omit(gene_list_interaction_skbr3)

# sort the list in decreasing order (#required for clusterProfiler)
gene_list_interaction_skbr3 = sort(gene_list_interaction_skbr3, decreasing = TRUE)

gse_interaction_skbr3 <- gseGO(geneList=gene_list_interaction_skbr3, 
                              ont ="ALL", 
                              keyType = "ENSEMBL", 
                              minGSSize = 3, 
                              maxGSSize = 800, 
                              pvalueCutoff = 0.05, 
                              verbose = TRUE, 
                              OrgDb = get('org.Hs.eg.db'), 
                              pAdjustMethod = "fdr", 
                              eps = 0,
                              nPermSimple = 10000)

require(DOSE)
dotplot(gse_interaction_skbr3, showCategory=20, split=".sign") + facet_grid(.~.sign)

set.seed(123)
x2 <- pairwise_termsim(gse_interaction_skbr3)
emapplot(x2)

p1 <- emapplot(x2, showCategory = 15, cex_label_category = 0.7)

p1 + theme_science(base_size = 8) + scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
ggsave("skbr3_interaction_emap.tiff", height = 5, width = 5, units = "in")


gse_interaction_skbr3 <- as.data.frame(gse_interaction_skbr3@result)

write.csv(gse_interaction_skbr3, "gse_interaction_skbr3.csv")

gse_interaction_skbr3 <- dplyr::arrange(gse_interaction_skbr3, desc(-p.adjust))

gse_interaction_skbr3_top  <- gse_interaction_skbr3[c(1:15),]

gse_interaction_skbr3_top$direction <- ifelse(gse_interaction_skbr3_top$NES > 0, "Enriched", "Depleted")

gse_int_skbr3 <- ggplot(gse_interaction_skbr3_top, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank(), panel.border = element_rect(fill = NA), 
        axis.line.y = element_line(colour= "black", size=0.5),
        axis.text.x = element_text(angle = 90, hjust = -1, vjust = 0.5))

gse_int_skbr3 + facet_wrap(~direction, ncol = 2, scales = "free_x")

ggsave("go_interaction_skbr3.tiff", height = 5, width = 5, units = "in")



######################################################
######
###### ADAR knockdown in skbr3
######
######
######################################################


adar_kd_skbr3 <- results(dds_skbr3, contrast=c("shrna1", "a", "s"))

adar_kd_skbr3_shrink <- lfcShrink(dds_skbr3, type="apeglm", coef = 2)

plotMA(adar_kd_skbr3)
plotMA(adar_kd_skbr3_shrink)

sum(adar_kd_skbr3$padj <0.01, na.rm=TRUE)
sum(adar_kd_skbr3_shrink$padj <0.01, na.rm=TRUE)

adar_kd_skbr3_df <- as.data.frame(adar_kd_skbr3)
adar_kd_skbr3_shrink_df <- as.data.frame(adar_kd_skbr3_shrink)

adar_kd_skbr3_df$id<- row.names(adar_kd_skbr3_df)
adar_kd_skbr3_shrink_df$id<- row.names(adar_kd_skbr3_shrink_df)

adar_kd_skbr3_df <- merge(adar_kd_skbr3_df, df_r, by = "id")
adar_kd_skbr3_shrink_df <- merge(adar_kd_skbr3_shrink_df, df_r, by = "id")

write.csv(adar_kd_skbr3_shrink_df, "adar_kd_skbr3_shrink_df.csv")

adar_kd_skbr3_shrink_df$Sig <- ifelse(adar_kd_skbr3_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(adar_kd_skbr3_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_adar_kd.tiff", height=4, width=4, unit="in")

adar_kd_skbr3_shrink_df$ISG <- grepl(paste(ISGS, collapse = "|"), adar_kd_skbr3_shrink_df$ensembl_gene_id)
adar_kd_skbr3_shrink_df$ATF4_targets <- grepl(paste(atf4_targets, collapse = "|"), adar_kd_skbr3_shrink_df$ensembl_gene_id)

ggplot(adar_kd_skbr3_shrink_df, aes(log2FoldChange, -log10(padj), colour = ISG, alpha = ISG)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_adar_kd_ISG.tiff", height=4, width=4, unit="in")


ggplot(adar_kd_skbr3_shrink_df, aes(log2FoldChange, -log10(padj), colour = ATF4_targets, alpha = ATF4_targets)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_adar_kd_ATF4.tiff", height=4, width=4, unit="in")


adar_kd_skbr3_shrink_df$rank <- -log10(adar_kd_skbr3_shrink_df$padj)*sign(adar_kd_skbr3_shrink_df$log2FoldChange)

adar_kd_skbr3_shrink_df <- do.call(data.frame, lapply(adar_kd_skbr3_shrink_df, function(x) replace(x, is.infinite(x), NA)))


# we want the log2 fold change 
gene_list_adar_kd_skbr3 <- adar_kd_skbr3_shrink_df$rank


# name the vector
names(gene_list_adar_kd_skbr3) <- adar_kd_skbr3_shrink_df$ensembl_gene_id

# omit any NA values 
gene_list_adar_kd_skbr3<-na.omit(gene_list_adar_kd_skbr3)

# sort the list in decreasing order (#required for clusterProfiler)
gene_list_adar_kd_skbr3 = sort(gene_list_adar_kd_skbr3, decreasing = TRUE)

gse_adar_kd_skbr3 <- gseGO(geneList=gene_list_adar_kd_skbr3, 
                          ont ="ALL", 
                          keyType = "ENSEMBL", 
                          minGSSize = 3, 
                          maxGSSize = 800, 
                          pvalueCutoff = 0.05, 
                          verbose = TRUE, 
                          OrgDb = get('org.Hs.eg.db'), 
                          pAdjustMethod = "fdr", 
                          eps = 0)

require(DOSE)
dotplot(gse_adar_kd_skbr3, showCategory=20, split=".sign") + facet_grid(.~.sign)

set.seed(123)
x2 <- pairwise_termsim(gse_adar_kd_skbr3)
emapplot(x2)

p1 <- emapplot(x2, showCategory = 15, cex_label_category = 0.7)

p1 + theme_science(base_size = 8) + scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
ggsave("skbr3_adar_kd_emap.tiff", height = 5, width = 5, units = "in")


gse_adar_kd_skbr3 <- as.data.frame(gse_adar_kd_skbr3@result)

write.csv(gse_adar_kd_skbr3, "gse_adar_kd_skbr3.csv")

gse_adar_kd_skbr3 <- dplyr::arrange(gse_adar_kd_skbr3, desc(-p.adjust))

gse_adar_kd_skbr3_top  <- gse_adar_kd_skbr3[c(1:15),]

gse_adar_kd_skbr3_top$direction <- ifelse(gse_adar_kd_skbr3_top$NES > 0, "Enriched", "Depleted")

gse_adar_skbr3 <- ggplot(gse_dhx9_kd_skbr3_top, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank(), panel.border = element_rect(fill = NA), axis.line.y = element_line(colour= "black", size=0.5))

gse_adar_skbr3 + facet_wrap(~direction, ncol = 2, scales = "free_x")

ggsave("go_adar_kd_skbr3.tiff", height = 5, width = 5, units = "in")

######################################################
######
###### DHX9 knockdown in skbr3
######
######
######################################################



dhx9_kd_skbr3 <- results(dds_skbr3, contrast=c("shrna2", "d", "s"))

dhx9_kd_skbr3_shrink <- lfcShrink(dds_skbr3, type="apeglm", coef = 3)

plotMA(dhx9_kd_skbr3)
plotMA(dhx9_kd_skbr3_shrink)

sum(dhx9_kd_skbr3$padj <0.01, na.rm=TRUE)
sum(dhx9_kd_skbr3_shrink$padj <0.01, na.rm=TRUE)

dhx9_kd_skbr3_df <- as.data.frame(dhx9_kd_skbr3)
dhx9_kd_skbr3_shrink_df <- as.data.frame(dhx9_kd_skbr3_shrink)

dhx9_kd_skbr3_df$id<- row.names(dhx9_kd_skbr3_df)
dhx9_kd_skbr3_shrink_df$id<- row.names(dhx9_kd_skbr3_shrink_df)

dhx9_kd_skbr3_df <- merge(dhx9_kd_skbr3_df, df_r, by = "id")
dhx9_kd_skbr3_shrink_df <- merge(dhx9_kd_skbr3_shrink_df, df_r, by = "id")

write.csv(dhx9_kd_skbr3_shrink_df, "dhx9_kd_skbr3_shrink_df.csv")

dhx9_kd_skbr3_shrink_df$Sig <- ifelse(dhx9_kd_skbr3_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(dhx9_kd_skbr3_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_dhx9_kd.tiff", height=4, width=4, unit="in")

dhx9_kd_skbr3_shrink_df$ISG <- grepl(paste(ISGS, collapse = "|"), dhx9_kd_skbr3_shrink_df$ensembl_gene_id)
dhx9_kd_skbr3_shrink_df$ATF4_targets <- grepl(paste(atf4_targets, collapse = "|"), dhx9_kd_skbr3_shrink_df$ensembl_gene_id)

ggplot(dhx9_kd_skbr3_shrink_df, aes(log2FoldChange, -log10(padj), colour = ISG, alpha = ISG)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_dhx9_kd_ISG.tiff", height=4, width=4, unit="in")


ggplot(dhx9_kd_skbr3_shrink_df, aes(log2FoldChange, -log10(padj), colour = ATF4_targets, alpha = ATF4_targets)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_dhx9_kd_ATF4.tiff", height=4, width=4, unit="in")


dhx9_kd_skbr3_shrink_df$rank <- -log10(dhx9_kd_skbr3_shrink_df$padj)*sign(dhx9_kd_skbr3_shrink_df$log2FoldChange)

dhx9_kd_skbr3_shrink_df <- do.call(data.frame, lapply(dhx9_kd_skbr3_shrink_df, function(x) replace(x, is.infinite(x), NA)))


# we want the log2 fold change 
gene_list_dhx9_kd_skbr3 <- dhx9_kd_skbr3_shrink_df$rank


# name the vector
names(gene_list_dhx9_kd_skbr3) <- dhx9_kd_skbr3_shrink_df$ensembl_gene_id

# omit any NA values 
gene_list_dhx9_kd_skbr3<-na.omit(gene_list_dhx9_kd_skbr3)

# sort the list in decreasing order (#required for clusterProfiler)
gene_list_dhx9_kd_skbr3 = sort(gene_list_dhx9_kd_skbr3, decreasing = TRUE)

gse_dhx9_kd_skbr3 <- gseGO(geneList=gene_list_dhx9_kd_skbr3, 
                          ont ="ALL", 
                          keyType = "ENSEMBL", 
                          minGSSize = 3, 
                          maxGSSize = 800, 
                          pvalueCutoff = 0.05, 
                          verbose = TRUE, 
                          OrgDb = get('org.Hs.eg.db'), 
                          pAdjustMethod = "fdr", 
                          eps = 0)

require(DOSE)
dotplot(gse_dhx9_kd_skbr3, showCategory=20, split=".sign") + facet_grid(.~.sign)

set.seed(123)
x2 <- pairwise_termsim(gse_dhx9_kd_skbr3)
emapplot(x2)

p1 <- emapplot(x2, showCategory = 15, cex_label_category = 0.7)

p1 + theme_science(base_size = 8) + scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
ggsave("skbr3_dhx9_kd_emap.tiff", height = 5, width = 5, units = "in")


gse_dhx9_kd_skbr3 <- as.data.frame(gse_dhx9_kd_skbr3@result)

write.csv(gse_dhx9_kd_skbr3, "gse_dhx9_kd_skbr3.csv")

gse_dhx9_kd_skbr3 <- dplyr::arrange(gse_dhx9_kd_skbr3, desc(-p.adjust))

gse_dhx9_kd_skbr3_top  <- gse_dhx9_kd_skbr3[c(1:15),]

gse_dhx9_kd_skbr3_top$direction <- ifelse(gse_dhx9_kd_skbr3_top$NES > 0, "Enriched", "Depleted")


gse_dhx9_skbr3 <- ggplot(gse_dhx9_kd_skbr3_top, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank(), panel.border = element_rect(fill = NA), axis.line.y = element_line(colour= "black", size=0.5))

gse_dhx9_skbr3 + facet_wrap(~direction, ncol = 2, scales = "free_x")

ggsave("go_dhx9_kd_skbr3.tiff", height = 5, width = 5, units = "in")


######################################################
######
###### INTERACTION in MCF7
######
######
######################################################


data_mcf7 <- as.matrix(df_r[,16:23])

row.names(data_mcf7) <- df_r$id

head(data_mcf7)

coldata_mcf7 <- DataFrame(shrna1 = factor(c("a", "a", "s", "s","a", "a", "s", "s")),
                          shrna2 = factor(c("d","s","d","s","d","s","d","s")),
                          row.names=as.character(colnames(data_mcf7)))
coldata_mcf7

ddsHTSeq_mcf7 <- DESeqDataSetFromMatrix(countData= data_mcf7, 
                                        colData = coldata_mcf7, 
                                        design = ~ shrna1 + shrna2 + shrna1:shrna2)
ddsHTSeq_mcf7

ddsHTSeq_mcf7$shrna1
ddsHTSeq_mcf7$shrna1 = relevel(ddsHTSeq_mcf7$shrna1, "s")
ddsHTSeq_mcf7$shrna1

ddsHTSeq_mcf7$shrna2
ddsHTSeq_mcf7$shrna2 = relevel(ddsHTSeq_mcf7$shrna2, "s")
ddsHTSeq_mcf7$shrna2

dds_mcf7 <- DESeq(ddsHTSeq_mcf7)

resultsNames(dds_mcf7)

keep <- rowSums(counts(dds_mcf7)) >= 10
dds_mcf7 <- dds_mcf7[keep,]


interaction_mcf7 <- results(dds_mcf7, name = "shrna1a.shrna2d")

interaction_mcf7_shrink <- lfcShrink(dds_mcf7, type="apeglm", coef = 4)

plotMA(interaction_mcf7)
plotMA(interaction_mcf7_shrink)

sum(interaction_mcf7$padj <0.01, na.rm=TRUE)
sum(interaction_mcf7_shrink$padj <0.01, na.rm=TRUE)

interaction_mcf7_df <- as.data.frame(interaction_mcf7)
interaction_mcf7_shrink_df <- as.data.frame(interaction_mcf7_shrink)

interaction_mcf7_df$id<- row.names(interaction_mcf7_df)
interaction_mcf7_shrink_df$id<- row.names(interaction_mcf7_shrink_df)

interaction_mcf7_df <- merge(interaction_mcf7_df, df_r, by = "id")
interaction_mcf7_shrink_df <- merge(interaction_mcf7_shrink_df, df_r, by = "id")

write.csv(interaction_mcf7_shrink_df, "interaction_mcf7_shrink_df.csv")


interaction_mcf7_shrink_df$Sig <- ifelse(interaction_mcf7_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(interaction_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_interaction.tiff", height=4, width=4, unit="in")

d <- plotCounts(dds_mcf7, gene=which.min(interaction_mcf7$padj), intgroup=c("shrna1", "shrna2"), returnData = TRUE)

base <- ggplot(d, aes(shrna1, log2(count+1), group = 1)) + 
  geom_point() + stat_summary(fun.y=mean, colour="red", geom="line") +
  xlab(NULL) + 
  ylab(NULL)

base + facet_wrap(~shrna2, ncol = 2)
ggsave("up_gene_ex.tiff")

d2 <- plotCounts(dds_mcf7, gene="34309_ENSG00000237423", intgroup=c("shrna1", "shrna2"), returnData = TRUE)

base <- ggplot(d2, aes(shrna1, log2(count+1), group = 1)) + 
  geom_point() + stat_summary(fun.y=mean, colour="red", geom="line") +
  xlab(NULL) + 
  ylab(NULL)

base + facet_wrap(~shrna2, ncol = 2)
ggsave("down_gene_ex.tiff")


d3 <- plotCounts(dds_mcf7, gene="7185_ENSG00000135829", intgroup=c("shrna1", "shrna2"), returnData = TRUE)

base <- ggplot(d3, aes(shrna1, log2(count+1), group = 1)) + 
  geom_point() + stat_summary(fun.y=mean, colour="red", geom="line") +
  xlab(NULL) + 
  ylab(NULL)

base + facet_wrap(~shrna2, ncol = 2)
ggsave("flat_gene_ex_dhx9.tiff")


d4 <- plotCounts(dds_mcf7, gene="10565_ENSG00000160710", intgroup=c("shrna1", "shrna2"), returnData = TRUE)

base <- ggplot(d4, aes(shrna1, log2(count+1), group = 1)) + 
  geom_point() + stat_summary(fun.y=mean, colour="red", geom="line") +
  xlab(NULL) + 
  ylab(NULL)

base + facet_wrap(~shrna2, ncol = 2)
ggsave("up_gene_ex_adar.tiff")




interaction_mcf7_shrink_df$ISG <- grepl(paste(ISGS, collapse = "|"), interaction_mcf7_shrink_df$ensembl_gene_id)
interaction_mcf7_shrink_df$ATF4_targets <- grepl(paste(atf4_targets, collapse = "|"), interaction_mcf7_shrink_df$ensembl_gene_id)
interaction_mcf7_shrink_df$oas <- grepl("ENSG00000111331|ENSG00000111335|ENSG00000089127", interaction_mcf7_shrink_df$ensembl_gene_id)


ggplot(interaction_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = ISG, alpha = ISG)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_interaction_ISG.tiff", height=4, width=4, unit="in")


ggplot(interaction_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = ATF4_targets, alpha = ATF4_targets)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_interaction_ATF4.tiff", height=4, width=4, unit="in")

ggplot(interaction_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = oas, alpha = oas)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_interaction_oas.tiff", height=4, width=4, unit="in")


interaction_mcf7_shrink_df$rank <- -log10(interaction_mcf7_shrink_df$padj)*sign(interaction_mcf7_shrink_df$log2FoldChange)

interaction_mcf7_shrink_df <- do.call(data.frame, lapply(interaction_mcf7_shrink_df, function(x) replace(x, is.infinite(x), NA)))


# we want the log2 fold change 
gene_list_interaction_mcf7 <- interaction_mcf7_shrink_df$rank

# name the vector
names(gene_list_interaction_mcf7) <- interaction_mcf7_shrink_df$ensembl_gene_id

# omit any NA values 
gene_list_interaction_mcf7<-na.omit(gene_list_interaction_mcf7)

# sort the list in decreasing order (#required for clusterProfiler)
gene_list_interaction_mcf7 = sort(gene_list_interaction_mcf7, decreasing = TRUE)

gse_interaction_mcf7 <- gseGO(geneList=gene_list_interaction_mcf7, 
                         ont ="ALL", 
                         keyType = "ENSEMBL", 
                         minGSSize = 3, 
                         maxGSSize = 800, 
                         pvalueCutoff = 0.05, 
                         verbose = TRUE, 
                         OrgDb = get('org.Hs.eg.db'), 
                         pAdjustMethod = "fdr", 
                         eps = 0,
                         nPermSimple = 10000)

require(DOSE)
dotplot(gse_interaction_mcf7, showCategory=20, split=".sign") + facet_grid(.~.sign)

set.seed(123)
x2 <- pairwise_termsim(gse_interaction_mcf7)
emapplot(x2)

p1 <- emapplot(x2, showCategory = 15, cex_label_category = 0.7)

p1 + theme_science(base_size = 8) + scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
ggsave("mcf7_interaction_emap.tiff", height = 5, width = 10, units = "in")

gse_interaction_mcf7 <- as.data.frame(gse_interaction_mcf7@result)

write.csv(gse_interaction_mcf7, "gse_interaction_mcf7.csv")

gse_interaction_mcf7 <- dplyr::arrange(gse_interaction_mcf7, desc(-p.adjust))

gse_interaction_mcf7_top  <- gse_interaction_mcf7[c(1:15),]

gse_interaction_mcf7_top$direction <- ifelse(gse_interaction_mcf7_top$NES > 0, "Enriched", "Depleted")

gse_int_mcf7 <- ggplot(gse_interaction_mcf7_top, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank(), panel.border = element_rect(fill = NA), 
        axis.line.y = element_line(colour= "black", size=0.5),
        axis.text.x = element_text(angle = 90, hjust = -1, vjust = 0.5))

gse_int_mcf7 + facet_wrap(~direction, ncol = 2, scales = "free_x")

ggsave("go_interaction_mcf7.tiff", height = 4, width = 7, units = "in")



######################################################
######
###### ADAR knockdown in MCF7
######
######
######################################################


adar_kd_mcf7 <- results(dds_mcf7, contrast=c("shrna1", "a", "s"))

adar_kd_mcf7_shrink <- lfcShrink(dds_mcf7, type="apeglm", coef = 2)

plotMA(adar_kd_mcf7)
plotMA(adar_kd_mcf7_shrink)

sum(adar_kd_mcf7$padj <0.01, na.rm=TRUE)
sum(adar_kd_mcf7_shrink$padj <0.01, na.rm=TRUE)

adar_kd_mcf7_df <- as.data.frame(adar_kd_mcf7)
adar_kd_mcf7_shrink_df <- as.data.frame(adar_kd_mcf7_shrink)

adar_kd_mcf7_df$id<- row.names(adar_kd_mcf7_df)
adar_kd_mcf7_shrink_df$id<- row.names(adar_kd_mcf7_shrink_df)

adar_kd_mcf7_df <- merge(adar_kd_mcf7_df, df_r, by = "id")
adar_kd_mcf7_shrink_df <- merge(adar_kd_mcf7_shrink_df, df_r, by = "id")

write.csv(adar_kd_mcf7_shrink_df, "adar_kd_mcf7_shrink_df.csv")

adar_kd_mcf7_shrink_df$Sig <- ifelse(adar_kd_mcf7_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(adar_kd_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_adar_kd.tiff", height=4, width=4, unit="in")

adar_kd_mcf7_shrink_df$ISG <- grepl(paste(ISGS, collapse = "|"), adar_kd_mcf7_shrink_df$ensembl_gene_id)
adar_kd_mcf7_shrink_df$ATF4_targets <- grepl(paste(atf4_targets, collapse = "|"), adar_kd_mcf7_shrink_df$ensembl_gene_id)

ggplot(adar_kd_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = ISG, alpha = ISG)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_adar_kd_ISG.tiff", height=4, width=4, unit="in")


ggplot(adar_kd_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = ATF4_targets, alpha = ATF4_targets)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_adar_kd_ATF4.tiff", height=4, width=4, unit="in")


adar_kd_mcf7_shrink_df$rank <- -log10(adar_kd_mcf7_shrink_df$padj)*sign(adar_kd_mcf7_shrink_df$log2FoldChange)

adar_kd_mcf7_shrink_df <- do.call(data.frame, lapply(adar_kd_mcf7_shrink_df, function(x) replace(x, is.infinite(x), NA)))


# we want the log2 fold change 
gene_list_adar_kd_mcf7 <- adar_kd_mcf7_shrink_df$rank


# name the vector
names(gene_list_adar_kd_mcf7) <- adar_kd_mcf7_shrink_df$ensembl_gene_id

# omit any NA values 
gene_list_adar_kd_mcf7<-na.omit(gene_list_adar_kd_mcf7)

# sort the list in decreasing order (#required for clusterProfiler)
gene_list_adar_kd_mcf7 = sort(gene_list_adar_kd_mcf7, decreasing = TRUE)

gse_adar_kd_mcf7 <- gseGO(geneList=gene_list_adar_kd_mcf7, 
                         ont ="ALL", 
                         keyType = "ENSEMBL", 
                         minGSSize = 3, 
                         maxGSSize = 800, 
                         pvalueCutoff = 0.05, 
                         verbose = TRUE, 
                         OrgDb = get('org.Hs.eg.db'), 
                         pAdjustMethod = "fdr", 
                         eps = 0)

require(DOSE)
dotplot(gse_adar_kd_mcf7, showCategory=20, split=".sign") + facet_grid(.~.sign)

set.seed(123)
x2 <- pairwise_termsim(gse_adar_kd_mcf7)
emapplot(x2)

p1 <- emapplot(x2, showCategory = 15, cex_label_category = 0.7)

p1 + theme_science(base_size = 8) + scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
ggsave("mcf7_adar_kd_emap.tiff", height = 5, width = 5, units = "in")


gse_adar_kd_mcf7 <- as.data.frame(gse_adar_kd_mcf7@result)

write.csv(gse_adar_kd_mcf7, "gse_adar_kd_mcf7.csv")

gse_adar_kd_mcf7 <- dplyr::arrange(gse_adar_kd_mcf7, desc(-p.adjust))

gse_adar_kd_mcf7_top  <- gse_adar_kd_mcf7[c(1:15),]

gse_adar_kd_mcf7_top$direction <- ifelse(gse_adar_kd_mcf7_top$NES > 0, "Enriched", "Depleted")

gse_adar_mcf7 <- ggplot(gse_adar_kd_mcf7_top, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank(), panel.border = element_rect(fill = NA), axis.line.y = element_line(colour= "black", size=0.5))

gse_adar_mcf7 + facet_wrap(~direction, ncol = 2, scales = "free_x")

ggsave("go_adar_kd_mcf7.tiff", height = 5, width = 5, units = "in")

######################################################
######
###### DHX9 knockdown in MCF7
######
######
######################################################



dhx9_kd_mcf7 <- results(dds_mcf7, contrast=c("shrna2", "d", "s"))

dhx9_kd_mcf7_shrink <- lfcShrink(dds_mcf7, type="apeglm", coef = 3)

plotMA(dhx9_kd_mcf7)
plotMA(dhx9_kd_mcf7_shrink)

sum(dhx9_kd_mcf7$padj <0.01, na.rm=TRUE)
sum(dhx9_kd_mcf7_shrink$padj <0.01, na.rm=TRUE)

dhx9_kd_mcf7_df <- as.data.frame(dhx9_kd_mcf7)
dhx9_kd_mcf7_shrink_df <- as.data.frame(dhx9_kd_mcf7_shrink)

dhx9_kd_mcf7_df$id<- row.names(dhx9_kd_mcf7_df)
dhx9_kd_mcf7_shrink_df$id<- row.names(dhx9_kd_mcf7_shrink_df)

dhx9_kd_mcf7_df <- merge(dhx9_kd_mcf7_df, df_r, by = "id")
dhx9_kd_mcf7_shrink_df <- merge(dhx9_kd_mcf7_shrink_df, df_r, by = "id")

write.csv(dhx9_kd_mcf7_shrink_df, "dhx9_kd_mcf7_shrink_df.csv")

dhx9_kd_mcf7_shrink_df$Sig <- ifelse(dhx9_kd_mcf7_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(dhx9_kd_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_dhx9_kd.tiff", height=4, width=4, unit="in")

dhx9_kd_mcf7_shrink_df$ISG <- grepl(paste(ISGS, collapse = "|"), dhx9_kd_mcf7_shrink_df$ensembl_gene_id)
dhx9_kd_mcf7_shrink_df$ATF4_targets <- grepl(paste(atf4_targets, collapse = "|"), dhx9_kd_mcf7_shrink_df$ensembl_gene_id)

ggplot(dhx9_kd_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = ISG, alpha = ISG)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_dhx9_kd_ISG.tiff", height=4, width=4, unit="in")


ggplot(dhx9_kd_mcf7_shrink_df, aes(log2FoldChange, -log10(padj), colour = ATF4_targets, alpha = ATF4_targets)) + geom_point() + scale_alpha_discrete(range = c(0.1,1)) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_dhx9_kd_ATF4.tiff", height=4, width=4, unit="in")


dhx9_kd_mcf7_shrink_df$rank <- -log10(dhx9_kd_mcf7_shrink_df$padj)*sign(dhx9_kd_mcf7_shrink_df$log2FoldChange)

dhx9_kd_mcf7_shrink_df <- do.call(data.frame, lapply(dhx9_kd_mcf7_shrink_df, function(x) replace(x, is.infinite(x), NA)))


# we want the log2 fold change 
gene_list_dhx9_kd_mcf7 <- dhx9_kd_mcf7_shrink_df$rank


# name the vector
names(gene_list_dhx9_kd_mcf7) <- dhx9_kd_mcf7_shrink_df$ensembl_gene_id

# omit any NA values 
gene_list_dhx9_kd_mcf7<-na.omit(gene_list_dhx9_kd_mcf7)

# sort the list in decreasing order (#required for clusterProfiler)
gene_list_dhx9_kd_mcf7 = sort(gene_list_dhx9_kd_mcf7, decreasing = TRUE)

gse_dhx9_kd_mcf7 <- gseGO(geneList=gene_list_dhx9_kd_mcf7, 
                         ont ="ALL", 
                         keyType = "ENSEMBL", 
                         minGSSize = 3, 
                         maxGSSize = 800, 
                         pvalueCutoff = 0.05, 
                         verbose = TRUE, 
                         OrgDb = get('org.Hs.eg.db'), 
                         pAdjustMethod = "fdr", 
                         eps = 0)

require(DOSE)
dotplot(gse_dhx9_kd_mcf7, showCategory=20, split=".sign") + facet_grid(.~.sign)

set.seed(123)
x2 <- pairwise_termsim(gse_dhx9_kd_mcf7)
emapplot(x2)

p1 <- emapplot(x2, showCategory = 15, cex_label_category = 0.7)

p1 + theme_science(base_size = 8) + scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
ggsave("mcf7_dhx9_kd_emap.tiff", height = 5, width = 5, units = "in")


gse_dhx9_kd_mcf7 <- as.data.frame(gse_dhx9_kd_mcf7@result)

write.csv(gse_dhx9_kd_mcf7, "gse_dhx9_kd_mcf7.csv")

gse_dhx9_kd_mcf7 <- dplyr::arrange(gse_dhx9_kd_mcf7, desc(-p.adjust))

gse_dhx9_kd_mcf7_top  <- gse_dhx9_kd_mcf7[c(1:15),]

gse_dhx9_kd_mcf7_top$direction <- ifelse(gse_dhx9_kd_mcf7_top$NES > 0, "Enriched", "Depleted")


gse_dhx9_mcf7 <- ggplot(gse_dhx9_kd_mcf7_top, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank(), panel.border = element_rect(fill = NA), axis.line.y = element_line(colour= "black", size=0.5))

gse_dhx9_mcf7 + facet_wrap(~direction, ncol = 2, scales = "free_x")

ggsave("go_dhx9_kd_mcf7.tiff", height = 5, width = 5, units = "in")




vst_mcf7 <- vst(dds_mcf7)

# Convert the DESeq transformed object to a data frame
vst_mcf7 <- assay(vst_mcf7)
vst_mcf7 <- as.data.frame(vst_mcf7)
colnames(vst_mcf7) <- gsub("sample.m", "MCF7_", colnames(vst_mcf7))


vst_skbr3 <- vst(dds_skbr3)

# Convert the DESeq transformed object to a data frame
vst_skbr3 <- assay(vst_skbr3)
vst_skbr3 <- as.data.frame(vst_skbr3)
colnames(vst_skbr3) <- gsub("sample.s", "SKBR3_", colnames(vst_skbr3))

vst_all <- merge(vst_mcf7, vst_skbr3, by ='row.names', all = TRUE)

row.names(vst_all) <- vst_all$Row.names

vst_all$Row.names <- NULL

vst_all <- as.data.frame(t(scale(t(vst_all))))

vst_all$id <- row.names(vst_all)

vst_all <- merge(vst_all, df_r_idents, by = "id")


row.names(vst_all) <- vst_all$id

vst_all$id <- NULL

vst_all <- na.omit(vst_all)


vst_all$ISG <- grepl(paste(ISGS, collapse = "|"), vst_all$ensembl_gene_id)

vst_all_ISG <- subset(vst_all, ISG == TRUE)

row.names(vst_all_ISG) <- vst_all_ISG$external_gene_name

vst_all_ISG$id <- NULL
vst_all_ISG$ensembl_gene_id <- NULL
vst_all_ISG$ISG <- NULL

vst_all_ISG <- na.omit(vst_all_ISG)

vst_all_ISG_long <- tidyr::gather(vst_all_ISG, sample, value, 1:16)

vst_all_ISG$external_gene_name <- NULL

vst_all_ISG_dendro <- as.dendrogram(hclust(d = dist(x = vst_all_ISG)))

# Create dendro
dendro_plot <- ggdendrogram(data = vst_all_ISG_dendro, rotate = TRUE)

# Preview the plot
print(dendro_plot)

vst_all_ISG_order <- order.dendrogram(vst_all_ISG_dendro)

vst_all_ISG$external_gene_name <- row.names(vst_all_ISG)

vst_all_ISG_long$external_gene_name <- factor(x = vst_all_ISG_long$external_gene_name,
                               levels = vst_all_ISG$external_gene_name[vst_all_ISG_order], 
                               ordered = TRUE)

vst_all_ISG_long$replicate <- str_sub(vst_all_ISG_long$sample, -1,-1)

vst_all_ISG_long$sample2 <- str_sub(vst_all_ISG_long$sample, end = -2)

vst_all_ISG_long$sample3 <- str_sub(vst_all_ISG_long$sample2,  start = -2)


vst_all_ISG_long$Cell_line <- ifelse(grepl("MCF7", vst_all_ISG_long$sample2), "MCF-7", "SK-BR-3")

# Make a heatmap
p1_isg <- ggplot(vst_all_ISG_long, aes(x=sample3, y=external_gene_name, fill=value, group = replicate)) + 
  geom_raster(position = "dodge") + 
  scale_fill_gradient2(low = palette[1], mid = 'white', high = palette[5]) + theme_science(base_size = 8) +
  theme(axis.text.y = element_text(size = 7), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_vline(xintercept = c(1.5,2.5,3.5,5.5,6.5,7.5), colour = "black", linetype = "dashed") +
  labs(x = "", fill = "RNA Expression\nz-Score") + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title.x = element_blank()) +
  scale_x_discrete(labels = c(
    "ad" = "shADAR\nshDHX9",
    "as" = "shADAR\nshSCR",
    "sd" = "shSCR\nshDHX9",
    "ss" = "shSCR\nshSCR"))
ggsave("isg_heatmap2.tiff", height = 8, width = 3, units = "in")

p1_isg + facet_wrap(~Cell_line, ncol = 2)

p1_isg + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Cell_line, nrow = 2, strip.position="left")

p2_isg <- ggplot(vst_all_ISG_long, aes(sample3, value)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA) + geom_quasirandom(alpha = 0.2) +  theme_science(base_size = 7) + 
  labs(y = "RNA Expression\n(z-score)", x = "") + 
  scale_x_discrete(labels = c(
    "ad" = "shADAR\nshDHX9",
    "as" = "shADAR\nshSCR",
    "sd" = "shSCR\nshDHX9",
    "ss" = "shSCR\nshSCR")) + geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
  theme(axis.text.y=element_blank(), 
        axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
        axis.title.x = element_text(vjust = -1,margin = unit(c(t = 0, r = -17, b = 0, l = 0), "mm")))
ggsave("isg_bp2.tiff", height= 5, width = 5, units = "in")

p2_isg

gA <- ggplotGrob(p2_isg + coord_flip() + facet_wrap(~Cell_line, nrow = 2) + theme(strip.text = element_blank()))
gB <- ggplotGrob(p1_isg + coord_flip() + 
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
                         legend.position = "top", strip.placement = "outside", strip.background = element_blank())+
                   facet_wrap(~Cell_line, nrow = 2, strip.position = "left"))

tiff("isg_heat_bp2_flip.tiff", height = 3.5, width = 6, units = "in", res = 300)
plot_grid(gB, gA, nrow = 1, align = "h", rel_widths = c(6,1), axis = "lrtb")
dev.off()






vst_all$atf4 <- grepl(paste(atf4_targets, collapse = "|"), vst_all$ensembl_gene_id)

vst_all_atf4 <- subset(vst_all, atf4 == TRUE)

row.names(vst_all_atf4) <- vst_all_atf4$external_gene_name

vst_all_atf4$id <- NULL
vst_all_atf4$ensembl_gene_id <- NULL
vst_all_atf4$atf4 <- NULL

vst_all_atf4 <- na.omit(vst_all_atf4)

vst_all_atf4_long <- tidyr::gather(vst_all_atf4, sample, value, 1:16)

vst_all_atf4$external_gene_name <- NULL

vst_all_atf4_dendro <- as.dendrogram(hclust(d = dist(x = vst_all_atf4)))

# Create dendro
dendro_plot <- ggdendrogram(data = vst_all_atf4_dendro, rotate = TRUE)

# Preview the plot
print(dendro_plot)

vst_all_atf4_order <- order.dendrogram(vst_all_atf4_dendro)

vst_all_atf4$external_gene_name <- row.names(vst_all_atf4)

vst_all_atf4_long$external_gene_name <- factor(x = vst_all_atf4_long$external_gene_name,
                                              levels = vst_all_atf4$external_gene_name[vst_all_atf4_order], 
                                              ordered = TRUE)

vst_all_atf4_long$replicate <- str_sub(vst_all_atf4_long$sample, -1,-1)

vst_all_atf4_long$sample2 <- str_sub(vst_all_atf4_long$sample, end = -2)

vst_all_atf4_long$sample3 <- str_sub(vst_all_atf4_long$sample2,  start = -2)

vst_all_atf4_long$Cell_line <- ifelse(grepl("MCF7", vst_all_atf4_long$sample2), "MCF-7", "SK-BR-3")

# Make a heatmap
p1_atf4 <- ggplot(vst_all_atf4_long, aes(x=sample3, y=external_gene_name, fill=value, group = replicate)) + 
  geom_raster(position = "dodge") + 
  scale_fill_gradient2(low = palette[1], mid = 'white', high = palette[5]) + theme_science(base_size = 8) +
  theme(axis.text.y = element_text(size = 7), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_vline(xintercept = c(1.5,2.5,3.5,5.5,6.5,7.5), colour = "black", linetype = "dashed") +
  labs(x = "", fill = "RNA Expression\nz-Score") + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title.x = element_blank()) +
  scale_x_discrete(labels = c(
    "ad" = "shADAR\nshDHX9",
    "as" = "shADAR\nshSCR",
    "sd" = "shSCR\nshDHX9",
    "ss" = "shSCR\nshSCR"))
ggsave("atf4_heatmap2.tiff", height = 8, width = 3, units = "in")

p1_atf4 + facet_wrap(~Cell_line, ncol = 2)

p1_atf4 + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Cell_line, nrow = 2, strip.position="left")

p2_atf4 <- ggplot(vst_all_atf4_long, aes(sample3, value)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA) + geom_quasirandom(alpha = 0.2) +  theme_science(base_size = 7) + 
  labs(y = "RNA Expression\n(z-score)", x = "") + 
  scale_x_discrete(labels = c(
    "ad" = "shADAR\nshDHX9",
    "as" = "shADAR\nshSCR",
    "sd" = "shSCR\nshDHX9",
    "ss" = "shSCR\nshSCR")) + geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
  theme(axis.text.y=element_blank(), 
        axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
        axis.title.x = element_text(vjust = -1,margin = unit(c(t = 0, r = -17, b = 0, l = 0), "mm")))
ggsave("atf4_bp2.tiff", height= 5, width = 5, units = "in")

p2_atf4

gA <- ggplotGrob(p2_atf4 + coord_flip() + facet_wrap(~Cell_line, nrow = 2) + theme(strip.text = element_blank()))
gB <- ggplotGrob(p1_atf4 + coord_flip() + 
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
                         legend.position = "top", strip.placement = "outside", strip.background = element_blank())+
                   facet_wrap(~Cell_line, nrow = 2, strip.position = "left"))

tiff("atf4_heat_bp2_flip.tiff", height = 3.5, width = 8, units = "in", res = 300)
plot_grid(gB, gA, nrow = 1, align = "h", rel_widths = c(6,1), axis = "lrtb")
dev.off()




vst_all$nfkb <- grepl(paste(nfkb, collapse = "|"), vst_all$ensembl_gene_id)

vst_all_nfkb <- subset(vst_all, nfkb == TRUE)

row.names(vst_all_nfkb) <- vst_all_nfkb$external_gene_name

vst_all_nfkb$id <- NULL
vst_all_nfkb$ensembl_gene_id <- NULL
vst_all_nfkb$nfkb <- NULL

vst_all_nfkb <- na.omit(vst_all_nfkb)

vst_all_nfkb_long <- tidyr::gather(vst_all_nfkb, sample, value, 1:16)

vst_all_nfkb$external_gene_name <- NULL

vst_all_nfkb_dendro <- as.dendrogram(hclust(d = dist(x = vst_all_nfkb)))

# Create dendro
dendro_plot <- ggdendrogram(data = vst_all_nfkb_dendro, rotate = TRUE)

# Preview the plot
print(dendro_plot)

vst_all_nfkb_order <- order.dendrogram(vst_all_nfkb_dendro)

vst_all_nfkb$external_gene_name <- row.names(vst_all_nfkb)

vst_all_nfkb_long$external_gene_name <- factor(x = vst_all_nfkb_long$external_gene_name,
                                               levels = vst_all_nfkb$external_gene_name[vst_all_nfkb_order], 
                                               ordered = TRUE)

vst_all_nfkb_long$replicate <- str_sub(vst_all_nfkb_long$sample, -1,-1)

vst_all_nfkb_long$sample2 <- str_sub(vst_all_nfkb_long$sample, end = -2)

vst_all_nfkb_long$sample3 <- str_sub(vst_all_nfkb_long$sample2,  start = -2)

vst_all_nfkb_long$Cell_line <- ifelse(grepl("MCF7", vst_all_nfkb_long$sample2), "MCF-7", "SK-BR-3")

# Make a heatmap
p1_nfkb <- ggplot(vst_all_nfkb_long, aes(x=sample3, y=external_gene_name, fill=value, group = replicate)) + 
  geom_raster(position = "dodge") + 
  scale_fill_gradient2(low = palette[1], mid = 'white', high = palette[5]) + theme_science(base_size = 8) +
  theme(axis.text.y = element_text(size = 7), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_vline(xintercept = c(1.5,2.5,3.5,5.5,6.5,7.5), colour = "black", linetype = "dashed") +
  labs(x = "", fill = "RNA Expression\nz-Score") + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title.x = element_blank()) +
  scale_x_discrete(labels = c(
    "ad" = "shADAR\nshDHX9",
    "as" = "shADAR\nshSCR",
    "sd" = "shSCR\nshDHX9",
    "ss" = "shSCR\nshSCR"))
ggsave("nfkb_heatmap2.tiff", height = 8, width = 3, units = "in")

p1_nfkb + facet_wrap(~Cell_line, ncol = 2)

p1_nfkb + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Cell_line, nrow = 2, strip.position="left")

p2_nfkb <- ggplot(vst_all_nfkb_long, aes(sample3, value)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA) + geom_quasirandom(alpha = 0.2) +  theme_science(base_size = 7) + 
  labs(y = "RNA Expression\n(z-score)", x = "") + 
  scale_x_discrete(labels = c(
    "ad" = "shADAR\nshDHX9",
    "as" = "shADAR\nshSCR",
    "sd" = "shSCR\nshDHX9",
    "ss" = "shSCR\nshSCR")) + geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
  theme(axis.text.y=element_blank(), 
        axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
        axis.title.x = element_text(vjust = -1,margin = unit(c(t = 0, r = -17, b = 0, l = 0), "mm")))
ggsave("nfkb_bp2.tiff", height= 5, width = 5, units = "in")

p2_nfkb

gA <- ggplotGrob(p2_nfkb + coord_flip() + facet_wrap(~Cell_line, nrow = 2) + theme(strip.text = element_blank()))
gB <- ggplotGrob(p1_nfkb + coord_flip() + 
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
                         legend.position = "top", strip.placement = "outside", strip.background = element_blank())+
                   facet_wrap(~Cell_line, nrow = 2, strip.position = "left"))

tiff("nfkb_heat_bp2_flip.tiff", height = 3.5, width = 13, units = "in", res = 300)
plot_grid(gB, gA, nrow = 1, align = "h", rel_widths = c(6,1), axis = "lrtb")
dev.off()




vst_all$upr <- grepl(paste(upr, collapse = "|"), vst_all$ensembl_gene_id)

vst_all_upr <- subset(vst_all, upr == TRUE)

row.names(vst_all_upr) <- vst_all_upr$external_gene_name

vst_all_upr$id <- NULL
vst_all_upr$ensembl_gene_id <- NULL
vst_all_upr$upr <- NULL

vst_all_upr <- na.omit(vst_all_upr)

vst_all_upr_long <- tidyr::gather(vst_all_upr, sample, value, 1:16)

vst_all_upr$external_gene_name <- NULL

vst_all_upr_dendro <- as.dendrogram(hclust(d = dist(x = vst_all_upr)))

# Create dendro
dendro_plot <- ggdendrogram(data = vst_all_upr_dendro, rotate = TRUE)

# Preview the plot
print(dendro_plot)

vst_all_upr_order <- order.dendrogram(vst_all_upr_dendro)

vst_all_upr$external_gene_name <- row.names(vst_all_upr)

vst_all_upr_long$external_gene_name <- factor(x = vst_all_upr_long$external_gene_name,
                                               levels = vst_all_upr$external_gene_name[vst_all_upr_order], 
                                               ordered = TRUE)

vst_all_upr_long$replicate <- str_sub(vst_all_upr_long$sample, -1,-1)

vst_all_upr_long$sample2 <- str_sub(vst_all_upr_long$sample, end = -2)


# Make a heatmap
p1_upr <- ggplot(vst_all_upr_long, aes(x=sample2, y=external_gene_name, fill=value, group = replicate)) + 
  geom_raster(position = "dodge") + 
  scale_fill_gradient2(low = palette[1], mid = 'white', high = palette[5]) + theme_science(base_size = 8) +
  theme(axis.text.y = element_text(size = 7), axis.title.x = element_blank()) +
  geom_vline(xintercept = 4.5, colour = "black") + 
  geom_vline(xintercept = c(1.5,2.5,3.5,5.5,6.5,7.5), colour = "black", linetype = "dashed") +
  labs(x = "SKBR3                               MCF7", fill = "RNA Expression\nz-Score") + 
  theme(legend.position = "bottom", plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_discrete(limits = rev, labels = c(
    "MCF7_ad" = "shADAR\nshDHX9",
    "MCF7_as" = "shADAR\nshSCR",
    "MCF7_sd" = "shSCR\nshDHX9",
    "MCF7_ss" = "shSCR\nshSCR",
    "SKBR3_ad" = "shADAR\nshDHX9",
    "SKBR3_as" = "shADAR\nshSCR",
    "SKBR3_sd" = "shSCR\nshDHX9",
    "SKBR3_ss" = "shSCR\nshSCR"))
ggsave("upr_heatmap2.tiff", height = 8, width = 3, units = "in")


p2_upr <- ggplot(vst_all_upr_long, aes(sample2, value)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA) + geom_quasirandom(alpha = 0.2) +  theme_science(base_size = 7) + 
  labs(y = "RNA Expression\n(z-score)", x = "") + geom_vline(xintercept = 4.5, colour = "black") +
  scale_x_discrete(limits = rev, labels = c(
    "MCF7_ad" = "shADAR\nshDHX9",
    "MCF7_as" = "shADAR\nshSCR",
    "MCF7_sd" = "shSCR\nshDHX9",
    "MCF7_ss" = "shSCR\nshSCR",
    "SKBR3_ad" = "shADAR\nshDHX9",
    "SKBR3_as" = "shADAR\nshSCR",
    "SKBR3_sd" = "shSCR\nshDHX9",
    "SKBR3_ss" = "shSCR\nshSCR")) + geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
  theme(axis.text.x=element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
        axis.title.y = element_text(vjust = -1,margin = unit(c(t = 0, r = -17, b = 0, l = 0), "mm")))
ggsave("upr_bp2.tiff", height= 5, width = 5, units = "in")



gA <- ggplotGrob(p2_upr)
gB <- ggplotGrob(p1_upr)

tiff("upr_heat_bp2.tiff", height = 12, width = 4.2, units = "in", res = 300)
plot_grid(gA, gB, ncol=1, align = "v", rel_heights = c(1,6))
dev.off()













#####All about those Alus

mad2 <- as.data.frame(fread("mad2.cntTable"))
row.names(mad2) <- mad2$`gene/TE`
mad2$`gene/TE` <- NULL
colnames(mad2) <- 'mad2'

mad1 <- as.data.frame(fread("mad1.cntTable"))
row.names(mad1) <- mad1$`gene/TE`
mad1$`gene/TE` <- NULL
colnames(mad1) <- 'mad1'

mas2 <- as.data.frame(fread("mas2.cntTable"))
row.names(mas2) <- mas2$`gene/TE`
mas2$`gene/TE` <- NULL
colnames(mas2) <- 'mas2'

mas1 <- as.data.frame(fread("mas1.cntTable"))
row.names(mas1) <- mas1$`gene/TE`
mas1$`gene/TE` <- NULL
colnames(mas1) <- 'mas1'


msd2 <- as.data.frame(fread("msd2.cntTable"))
row.names(msd2) <- msd2$`gene/TE`
msd2$`gene/TE` <- NULL
colnames(msd2) <- 'msd2'

msd1 <- as.data.frame(fread("msd1.cntTable"))
row.names(msd1) <- msd1$`gene/TE`
msd1$`gene/TE` <- NULL
colnames(msd1) <- 'msd1'


mss2 <- as.data.frame(fread("mss2.cntTable"))
row.names(mss2) <- mss2$`gene/TE`
mss2$`gene/TE` <- NULL
colnames(mss2) <- 'mss2'

mss1 <- as.data.frame(fread("mss1.cntTable"))
row.names(mss1) <- mss1$`gene/TE`
mss1$`gene/TE` <- NULL
colnames(mss1) <- 'mss1'




sad2 <- as.data.frame(fread("sad2.cntTable"))
row.names(sad2) <- sad2$`gene/TE`
sad2$`gene/TE` <- NULL
colnames(sad2) <- 'sad2'

sad1 <- as.data.frame(fread("sad1.cntTable"))
row.names(sad1) <- sad1$`gene/TE`
sad1$`gene/TE` <- NULL
colnames(sad1) <- 'sad1'

sas2 <- as.data.frame(fread("sas2.cntTable"))
row.names(sas2) <- sas2$`gene/TE`
sas2$`gene/TE` <- NULL
colnames(sas2) <- 'sas2'

sas1 <- as.data.frame(fread("sas1.cntTable"))
row.names(sas1) <- sas1$`gene/TE`
sas1$`gene/TE` <- NULL
colnames(sas1) <- 'sas1'


ssd2 <- as.data.frame(fread("ssd2.cntTable"))
row.names(ssd2) <- ssd2$`gene/TE`
ssd2$`gene/TE` <- NULL
colnames(ssd2) <- 'ssd2'

ssd1 <- as.data.frame(fread("ssd1.cntTable"))
row.names(ssd1) <- ssd1$`gene/TE`
ssd1$`gene/TE` <- NULL
colnames(ssd1) <- 'ssd1'


sss2 <- as.data.frame(fread("sss2.cntTable"))
row.names(sss2) <- sss2$`gene/TE`
sss2$`gene/TE` <- NULL
colnames(sss2) <- 'sss2'

sss1 <- as.data.frame(fread("sss1.cntTable"))
row.names(sss1) <- sss1$`gene/TE`
sss1$`gene/TE` <- NULL
colnames(sss1) <- 'sss1'




mcf7_te <- cbind(mad2,mad1,mas2,mas1, msd2, msd1, mss2, mss1)

skbr3_te <- cbind(sad2,sad1,sas2,sas1, ssd2, ssd1, sss2, sss1)

all_te <- cbind(mcf7_te, skbr3_te)

colnames(all_te) <- paste("sample.", colnames(all_te), sep = "")

#write.table(all_te, file='all_te.tsv', quote=FALSE, sep='\t')


data_mcf7_te <- as.matrix(mcf7_te)

head(data_mcf7_te)

coldata_mcf7_te <- DataFrame(shrna1 = factor(c("a", "a", "a", "a","s", "s", "s", "s")),
                           shrna2 = factor(c("d","d","s","s","d","d","s","s")),
                           row.names=as.character(colnames(data_mcf7_te)))
coldata_mcf7_te

ddsHTSeq_mcf7_te <- DESeqDataSetFromMatrix(countData= data_mcf7_te, 
                                         colData = coldata_mcf7_te, 
                                         design = ~ shrna1 + shrna2 + shrna1:shrna2 )
ddsHTSeq_mcf7_te

ddsHTSeq_mcf7_te$shrna1
ddsHTSeq_mcf7_te$shrna1 = relevel(ddsHTSeq_mcf7_te$shrna1, "s")
ddsHTSeq_mcf7_te$shrna1

ddsHTSeq_mcf7_te$shrna2
ddsHTSeq_mcf7_te$shrna2 = relevel(ddsHTSeq_mcf7_te$shrna2, "s")
ddsHTSeq_mcf7_te$shrna2

dds_mcf7_te <- DESeq(ddsHTSeq_mcf7_te)

keep <- rowSums(counts(dds_mcf7_te)) >= 10
dds_mcf7_te <- dds_mcf7_te[keep,]

resultsNames(dds_mcf7_te)

interactionmcf7_te <- results(dds_mcf7_te, name = "shrna1a.shrna2d")

interactionmcf7_te_shrink <- lfcShrink(dds_mcf7_te, type="apeglm", coef = 4)

plotMA(interactionmcf7_te)
plotMA(interactionmcf7_te_shrink)

sum(interactionmcf7_te$padj <0.01, na.rm=TRUE)
sum(interactionmcf7_te_shrink$padj <0.01, na.rm=TRUE)

interactionmcf7_te_df <- as.data.frame(interactionmcf7_te)
interactionmcf7_te_shrink_df <- as.data.frame(interactionmcf7_te_shrink)

interactionmcf7_te_df$id<- row.names(interactionmcf7_te_df)
interactionmcf7_te_shrink_df$id<- row.names(interactionmcf7_te_shrink_df)

interactionmcf7_te_shrink_df$Sig <- ifelse(interactionmcf7_te_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(interactionmcf7_te_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 

te_ident <- function(x) { 
  if(!grepl("SINE|LINE|Retro|LTR", x)) y <- "Absent"
  if(grepl("SINE", x)) y <- "SINE/Alu"
  if(grepl("LINE", x)) y <- "LINE"
  if(grepl("LTR", x)) y <- "LTR"
  if(grepl("Retro", x)) y <- "Retrotransposon"
  return(y)
}

#applies Amber_gln function
interactionmcf7_te_shrink_df$TE <- sapply(interactionmcf7_te_shrink_df$id, te_ident)

ggplot(subset(interactionmcf7_te_shrink_df, !TE == "Absent"), aes(log2FoldChange, -log10(padj), colour = TE)) + geom_point(alpha = 0.3) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = c(palette[1], palette[3], palette[5], palette[7])) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_interaction_te.tiff", height=4, width=5, unit="in")

interactionmcf7_te_shrink_df_rm <- interactionmcf7_te_shrink_df[!grepl("ENSG",interactionmcf7_te_shrink_df$id),]

write.csv(interactionmcf7_te_shrink_df_rm, "interaction_mcf7_te.csv")


dhx9_kd_mcf7_te <- results(dds_mcf7_te, contrast = c("shrna2", "d", "s"))

dhx9_kd_mcf7_te_shrink <- lfcShrink(dds_mcf7_te, type="apeglm", coef = 3)

plotMA(dhx9_kd_mcf7_te)
plotMA(dhx9_kd_mcf7_te_shrink)

sum(dhx9_kd_mcf7_te$padj <0.01, na.rm=TRUE)
sum(dhx9_kd_mcf7_te_shrink$padj <0.01, na.rm=TRUE)

dhx9_kd_mcf7_te_df <- as.data.frame(dhx9_kd_mcf7_te)
dhx9_kd_mcf7_te_shrink_df <- as.data.frame(dhx9_kd_mcf7_te_shrink)

dhx9_kd_mcf7_te_df$id<- row.names(dhx9_kd_mcf7_te_df)
dhx9_kd_mcf7_te_shrink_df$id<- row.names(dhx9_kd_mcf7_te_shrink_df)

dhx9_kd_mcf7_te_shrink_df$Sig <- ifelse(dhx9_kd_mcf7_te_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(dhx9_kd_mcf7_te_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 

#applies function to id class of TE
dhx9_kd_mcf7_te_shrink_df$TE <- sapply(dhx9_kd_mcf7_te_shrink_df$id, te_ident)

ggplot(subset(dhx9_kd_mcf7_te_shrink_df, !TE == "Absent"), aes(log2FoldChange, -log10(padj), colour = TE)) + geom_point(alpha = 0.3) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = c(palette[1], palette[3], palette[5], palette[7])) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_dhx9_kd_te.tiff", height=4, width=5, unit="in")

dhx9_kd_mcf7_te_shrink_df_rm <- dhx9_kd_mcf7_te_shrink_df[!grepl("ENSG",dhx9_kd_mcf7_te_shrink_df$id),]

write.csv(dhx9_kd_mcf7_te_shrink_df_rm, "dhx9_kd_mcf7_te.csv")


adar_kd_mcf7_te <- results(dds_mcf7_te, contrast = c("shrna1", "a", "s"))

adar_kd_mcf7_te_shrink <- lfcShrink(dds_mcf7_te, type="apeglm", coef = 2)

plotMA(adar_kd_mcf7_te)
plotMA(adar_kd_mcf7_te_shrink)

sum(adar_kd_mcf7_te$padj <0.01, na.rm=TRUE)
sum(adar_kd_mcf7_te_shrink$padj <0.01, na.rm=TRUE)

adar_kd_mcf7_te_df <- as.data.frame(adar_kd_mcf7_te)
adar_kd_mcf7_te_shrink_df <- as.data.frame(adar_kd_mcf7_te_shrink)

adar_kd_mcf7_te_df$id<- row.names(adar_kd_mcf7_te_df)
adar_kd_mcf7_te_shrink_df$id<- row.names(adar_kd_mcf7_te_shrink_df)

adar_kd_mcf7_te_shrink_df$Sig <- ifelse(adar_kd_mcf7_te_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(adar_kd_mcf7_te_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 

#applies function to id class of TE
adar_kd_mcf7_te_shrink_df$TE <- sapply(adar_kd_mcf7_te_shrink_df$id, te_ident)

ggplot(subset(adar_kd_mcf7_te_shrink_df, !TE == "Absent"), aes(log2FoldChange, -log10(padj), colour = TE)) + geom_point(alpha = 0.3) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = c(palette[1], palette[3], palette[5], palette[7])) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_mcf7_adar_kd_te.tiff", height=4, width=5, unit="in")


adar_kd_mcf7_te_shrink_df_rm <- adar_kd_mcf7_te_shrink_df[!grepl("ENSG",adar_kd_mcf7_te_shrink_df$id),]

write.csv(adar_kd_mcf7_te_shrink_df_rm, "adar_kd_mcf7_te.csv")



data_skbr3_te <- as.matrix(skbr3_te)

head(data_skbr3_te)

coldata_skbr3_te <- DataFrame(shrna1 = factor(c("a", "a", "a", "a","s", "s", "s", "s")),
                             shrna2 = factor(c("d","d","s","s","d","d","s","s")),
                             row.names=as.character(colnames(data_skbr3_te)))
coldata_skbr3_te

ddsHTSeq_skbr3_te <- DESeqDataSetFromMatrix(countData= data_skbr3_te, 
                                           colData = coldata_skbr3_te, 
                                           design = ~ shrna1 + shrna2 + shrna1:shrna2 )
ddsHTSeq_skbr3_te

ddsHTSeq_skbr3_te$shrna1
ddsHTSeq_skbr3_te$shrna1 = relevel(ddsHTSeq_skbr3_te$shrna1, "s")
ddsHTSeq_skbr3_te$shrna1

ddsHTSeq_skbr3_te$shrna2
ddsHTSeq_skbr3_te$shrna2 = relevel(ddsHTSeq_skbr3_te$shrna2, "s")
ddsHTSeq_skbr3_te$shrna2

dds_skbr3_te <- DESeq(ddsHTSeq_skbr3_te)

keep <- rowSums(counts(dds_skbr3_te)) >= 10
dds_skbr3_te <- dds_skbr3_te[keep,]

resultsNames(dds_skbr3_te)

interactionskbr3_te <- results(dds_skbr3_te, name = "shrna1a.shrna2d")

interactionskbr3_te_shrink <- lfcShrink(dds_skbr3_te, type="apeglm", coef = 4)

plotMA(interactionskbr3_te)
plotMA(interactionskbr3_te_shrink)

sum(interactionskbr3_te$padj <0.01, na.rm=TRUE)
sum(interactionskbr3_te_shrink$padj <0.01, na.rm=TRUE)

interactionskbr3_te_df <- as.data.frame(interactionskbr3_te)
interactionskbr3_te_shrink_df <- as.data.frame(interactionskbr3_te_shrink)

interactionskbr3_te_df$id<- row.names(interactionskbr3_te_df)
interactionskbr3_te_shrink_df$id<- row.names(interactionskbr3_te_shrink_df)

interactionskbr3_te_shrink_df$Sig <- ifelse(interactionskbr3_te_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(interactionskbr3_te_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 

te_ident <- function(x) { 
  if(!grepl("SINE|LINE|Retro|LTR", x)) y <- "Absent"
  if(grepl("SINE", x)) y <- "SINE/Alu"
  if(grepl("LINE", x)) y <- "LINE"
  if(grepl("LTR", x)) y <- "LTR"
  if(grepl("Retro", x)) y <- "Retrotransposon"
  return(y)
}

#applies function to id class of TE
interactionskbr3_te_shrink_df$TE <- sapply(interactionskbr3_te_shrink_df$id, te_ident)

ggplot(subset(interactionskbr3_te_shrink_df, !TE == "Absent"), aes(log2FoldChange, -log10(padj), colour = TE)) + geom_point(alpha = 0.3) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = c(palette[1], palette[3], palette[5], palette[7])) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_interaction_te.tiff", height=4, width=5, unit="in")

interactionskbr3_te_shrink_df_rm <- interactionskbr3_te_shrink_df[!grepl("ENSG",interactionskbr3_te_shrink_df$id),]

write.csv(interactionskbr3_te_shrink_df_rm, "interaction_skbr3_te.csv")


dhx9_kd_skbr3_te <- results(dds_skbr3_te, contrast = c("shrna2", "d", "s"))

dhx9_kd_skbr3_te_shrink <- lfcShrink(dds_skbr3_te, type="apeglm", coef = 3)

plotMA(dhx9_kd_skbr3_te)
plotMA(dhx9_kd_skbr3_te_shrink)

sum(dhx9_kd_skbr3_te$padj <0.01, na.rm=TRUE)
sum(dhx9_kd_skbr3_te_shrink$padj <0.01, na.rm=TRUE)

dhx9_kd_skbr3_te_df <- as.data.frame(dhx9_kd_skbr3_te)
dhx9_kd_skbr3_te_shrink_df <- as.data.frame(dhx9_kd_skbr3_te_shrink)

dhx9_kd_skbr3_te_df$id<- row.names(dhx9_kd_skbr3_te_df)
dhx9_kd_skbr3_te_shrink_df$id<- row.names(dhx9_kd_skbr3_te_shrink_df)

dhx9_kd_skbr3_te_shrink_df$Sig <- ifelse(dhx9_kd_skbr3_te_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(dhx9_kd_skbr3_te_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 

#applies function to id class of TE
dhx9_kd_skbr3_te_shrink_df$TE <- sapply(dhx9_kd_skbr3_te_shrink_df$id, te_ident)

ggplot(subset(dhx9_kd_skbr3_te_shrink_df, !TE == "Absent"), aes(log2FoldChange, -log10(padj), colour = TE)) + geom_point(alpha = 0.3) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = c(palette[1], palette[3], palette[5], palette[7])) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_dhx9_kd_te.tiff", height=4, width=5, unit="in")

dhx9_kd_skbr3_te_shrink_df_rm <- dhx9_kd_skbr3_te_shrink_df[!grepl("ENSG",dhx9_kd_skbr3_te_shrink_df$id),]

write.csv(dhx9_kd_skbr3_te_shrink_df_rm, "dhx9_kd_skbr3_te.csv")



adar_kd_skbr3_te <- results(dds_skbr3_te, contrast = c("shrna1", "a", "s"))

adar_kd_skbr3_te_shrink <- lfcShrink(dds_skbr3_te, type="apeglm", coef = 2)

plotMA(adar_kd_skbr3_te)
plotMA(adar_kd_skbr3_te_shrink)

sum(adar_kd_skbr3_te$padj <0.01, na.rm=TRUE)
sum(adar_kd_skbr3_te_shrink$padj <0.01, na.rm=TRUE)

adar_kd_skbr3_te_df <- as.data.frame(adar_kd_skbr3_te)
adar_kd_skbr3_te_shrink_df <- as.data.frame(adar_kd_skbr3_te_shrink)

adar_kd_skbr3_te_df$id<- row.names(adar_kd_skbr3_te_df)
adar_kd_skbr3_te_shrink_df$id<- row.names(adar_kd_skbr3_te_shrink_df)

adar_kd_skbr3_te_shrink_df$Sig <- ifelse(adar_kd_skbr3_te_shrink_df$padj<0.05, "<0.05", ">0.05")

ggplot(adar_kd_skbr3_te_shrink_df, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = palette2) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 

#applies function to id class of TE
adar_kd_skbr3_te_shrink_df$TE <- sapply(adar_kd_skbr3_te_shrink_df$id, te_ident)

ggplot(subset(adar_kd_skbr3_te_shrink_df, !TE == "Absent"), aes(log2FoldChange, -log10(padj), colour = TE)) + geom_point(alpha = 0.3) +
  theme_science()  + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + scale_colour_manual(values = c(palette[1], palette[3], palette[5], palette[7])) +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey") + scale_x_continuous(limits = c(-8,8), breaks = c(-5,-2.5,0,2.5,5)) 
ggsave("volcano_skbr3_adar_kd_te.tiff", height=4, width=5, unit="in")

adar_kd_skbr3_te_shrink_df_rm <- adar_kd_skbr3_te_shrink_df[!grepl("ENSG",adar_kd_skbr3_te_shrink_df$id),]

write.csv(adar_kd_skbr3_te_shrink_df_rm, "adar_kd_skbr3_te.csv")










