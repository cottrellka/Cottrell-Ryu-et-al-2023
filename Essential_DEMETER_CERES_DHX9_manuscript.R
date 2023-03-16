library(ggplot2)
library(edgeR)
library(dplyr)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(data.table)
library(limma)
library(Biobase)
library(DESeq2)
library(biomaRt)
library(ggpubr)
library(MetBrewer)

font_import()

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

#this function stacks ggplot graphs, downloaded from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#used with scale_axis_continuous to only show two decimal points, downloaded from https://code-examples.net/en/q/24eda9a
scaleFUN <- function(x) sprintf("%.2f", x)

#color blind friendly palette
palette <- met.brewer(name="Veronese", n=7, type="discrete")

met.brewer(name="Veronese", n=7, type="discrete")

palette2 <- c(palette[7], palette[3])
palette3 <- c(palette[7], palette[5], palette[3])

setwd("/Users/cottr/Box Sync/BRCA/Essential Genes/")

#cell line annotations from Marcotte et al., 2016, https://github.com/neellab/bfg/tree/gh-pages/data/annotations
cell_line_subtypes <- read.delim("cell_line_subtypes_corrected.txt")

#read in DEMETER2 scores from DepMap portal, get just ADAR scores
depmap <- fread("D2_combined_gene_dep_scores.csv")
depmap <- subset(depmap, V1 == "DHX9 (1660)")
#make row.names = ADAR (103)
depmap <- data.frame(depmap, row.names = depmap$V1)
#get rid of column 1
depmap$V1 <- NULL
#transpose, rows are now cells and the one column is ADAR DEMETER2
depmap <- as.data.frame(t(depmap))
#change column name
colnames(depmap) <- "DEMETER2_DHX9"

#split row.names by '_' and make new columns with cell_line and Site
out <- strsplit(as.character(row.names(depmap)), "_", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

depmap$cell_line <- out$V1
depmap$Site <- out$V2


#Choose only Breast cell lines and get rid of Site
depmap <- subset(depmap, Site == "BREAST")
depmap$Site <- NULL

#read in CHRONOS data and annotations for mapping
CHRONOS <- fread("CRISPR_gene_effect.csv")
Cell_annotations <- fread("sample_info.csv")

#New data.frame for CHRONOS with just ADAR scores and cell lines, merge with annotations
CHRONOS <- data.frame(CHRONOS = CHRONOS$`DHX9 (1660)`, DepMap_ID = CHRONOS$DepMap_ID)

CHRONOS <- merge(CHRONOS, Cell_annotations, by = "DepMap_ID", all = TRUE)

#Choose only Breast cell lines
CHRONOS <- subset(CHRONOS, lineage == "breast")

#merge demeter2 and CHRONOS scores
Dependency <- merge(depmap, CHRONOS, by.x = "cell_line", by.y = "stripped_cell_line_name", all = TRUE)

Dependency$cell_line <- tolower(Dependency$cell_line)

#merge with cell line annotations
Dependency <- merge(cell_line_subtypes, Dependency, by = "cell_line", all = TRUE)

#Make plots of interest
ggplot(subset(Dependency, !is.na(DEMETER2_DHX9) & !is.na(subtype_three_receptor)), aes(reorder(cell_line, -DEMETER2_DHX9), DEMETER2_DHX9, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_x_discrete(position = "top") +
  scale_fill_manual(values = palette3) + theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="DHX9 DEMETER2 Score", fill = "Subtype")
ggsave("DHX9_BCCL_Depmap.tiff", height = 3, width =5, units = "in")

ggplot(subset(Dependency, !is.na(CHRONOS) & !is.na(subtype_three_receptor)), aes(reorder(cell_line, -CHRONOS), CHRONOS, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_x_discrete(position = "top") +
  scale_fill_manual(values = palette3) + theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="DHX9 CHRONOS Score", fill = "Subtype")
ggsave("DHX9_BCCL_CHRONOS.tiff")




#read in CCLE rnaseq data, values are log2(TPM+1)
CCLE <- fread("CCLE_expression_full.csv")

CCLE <- merge(CCLE, Cell_annotations, by.x = "V1", by.y = "DepMap_ID")

CCLE_breast <- subset(CCLE, lineage == "breast")

CCLE_breast <- CCLE_breast[,1:53971]



#merge with dependency data
Dependency_CCLE <- merge(Dependency, CCLE_breast, by.x = "DepMap_ID", by.y = "V1")

#make row.names cell_line names
row.names(Dependency_CCLE) <- Dependency_CCLE$cell_line

#make plots of interest

ggplot(CCLE, aes(`ADAR (ENSG00000160710)`, `DHX9 (ENSG00000135829)`, colour = lineage)) + 
  geom_point() +
  theme_science() +
  labs(x = "ADAR", y="DHX9") + geom_smooth(method = "lm")
ggsave("DHX9_ADAR_all_lineage.tiff", height = 6, width = 15, units = "in")
cor.test(CCLE$`ADAR (ENSG00000160710)`, CCLE$`DHX9 (ENSG00000135829)`)

ggplot(CCLE, aes(`ADAR (ENSG00000160710)`, `DHX9 (ENSG00000135829)`)) + 
  geom_point() +
  theme_science() +
  labs(x = "ADAR", y="DHX9") + geom_smooth(method = "lm")
ggsave("DHX9_ADAR_all.tiff", height = 3.7, width = 3.5, units = "in")
cor.test(CCLE$`ADAR (ENSG00000160710)`, CCLE$`DHX9 (ENSG00000135829)`)

ggplot(Dependency_CCLE, aes(`ADAR (ENSG00000160710)`, `DHX9 (ENSG00000135829)`)) + 
  geom_point() + theme_science() +
  geom_smooth(method = "lm", colour = palette2[2]) + labs(y = "ADAR Expression\n (log2(TPM+1))", x = "DHX9 Expression\n (log2(TPM+1))", colour = "") +
  theme(legend.position = "bottom") + stat_cor(method="pearson", alpha = 1)
ggsave("DHX9_ADAR_breast.tiff", height = 3.7, width = 3.5, units = "in")
cor.test(Dependency_CCLE$`ADAR (ENSG00000160710)`, Dependency_CCLE$`DHX9 (ENSG00000135829)`)


ggplot(Dependency_CCLE, aes(`ADAR (ENSG00000160710)`, `DDX17 (ENSG00000100201)`)) + 
  geom_point() +
  theme_science() +
  labs(x = "ADAR", y="DDX17") + geom_smooth(method = "lm")
cor.test(Dependency_CCLE$`ADAR (ENSG00000160710)`, Dependency_CCLE$`DDX17 (ENSG00000100201)`)

ggplot(Dependency_CCLE, aes(`ADAR (ENSG00000160710)`, `DDX54 (ENSG00000123064)`)) + 
  geom_point() +
  theme_science() +
  labs(x = "ADAR", y="DDX17") + geom_smooth(method = "lm")
cor.test(Dependency_CCLE$`ADAR (ENSG00000160710)`, Dependency_CCLE$`DDX54 (ENSG00000123064)`)


#find correlation between ADAR expression and expression of all other genes, pearson
CCLE_breast <- as.data.frame(CCLE_breast)

CORS_rna_p <- matrix(nrow=53996, ncol=3) 
for (i in 2:length(CCLE_breast)) {
  a <- cor.test(CCLE_breast$`ADAR (ENSG00000160710)`, CCLE_breast[,i], method = "pearson")
  CORS_rna_p[i,] <- c(colnames(CCLE_breast)[i], a$estimate, a$p.value)
}

CORS_rna_p <- as.data.frame(na.omit(CORS_rna_p))

CORS_rna_p$CC <- as.numeric(as.character(CORS_rna_p$V2))
CORS_rna_p$p_value <- as.numeric(as.character(CORS_rna_p$V3))
CORS_rna_p$FDR <- p.adjust(CORS_rna_p$p_value, method = "fdr", n = length(CORS_rna_p$p_value))
CORS_rna_p$method = "Pearson"

#find correlation between ADAR expression and expression of all other genes, spearman

CORS_rna_s <- matrix(nrow=53996, ncol=3) 
for (i in 2:length(CCLE_breast)) {
  a <- cor.test(CCLE_breast$`ADAR (ENSG00000160710)`, CCLE_breast[,i], method = "spearman")
  CORS_rna_s[i,] <- c(colnames(CCLE_breast)[i], a$estimate, a$p.value)
}

CORS_rna_s <- as.data.frame(na.omit(CORS_rna_s))

CORS_rna_s$CC <- as.numeric(as.character(CORS_rna_s$V2))
CORS_rna_s$p_value <- as.numeric(as.character(CORS_rna_s$V3))
CORS_rna_s$FDR <- p.adjust(CORS_rna_s$p_value, method = "fdr", n = length(CORS_rna_s$p_value))
CORS_rna_s$method = "Spearman"

CORS_rna_merge <- rbind(CORS_rna_p, CORS_rna_s)

#subset for only DDX## and DHX## genes

CORS_rna_helicase <- subset(CORS_rna_merge, grepl("DDX|DHX", CORS_rna_merge$V1))

CORS_rna_helicase$rank <- CORS_rna_helicase$Pearson + CORS_rna_helicase$Spearman /2

out <- strsplit(CORS_rna_helicase$V1, " ", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

CORS_rna_helicase$Gene <- out$V1

#subset for only helicases identified by APEX labeling

CORS_rna_helicase <- subset(CORS_rna_helicase, Gene %in% c("DHX9","DDX17", "DDX24", "DDX18", "DDX21", "DDX5", "DDX54", "DDX20", "DDX27"))

ggplot(CORS_rna_helicase, aes(CC, reorder(Gene, CC), size = -log10(FDR), colour = method)) + 
  geom_point() + theme_science() + labs(x = "Correlation Coefficient", y = "", colour = "") +
  scale_colour_manual(values = palette2)
ggsave("cors_ccle.tiff", height = 3.5, width = 4.5)

write.csv(CORS_rna_helicase, "CORS_rna_helicase.csv")




                            
#read CCLE rnaseq transcript data - pre-normalized
CCLE_transcripts <- fread("CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt")

#get only ADAR transcripts
CCLE_transcripts_ADAR <- subset(CCLE_transcripts, grepl("ENSG00000160710", CCLE_transcripts$gene_id))

#make row.names = transcript_id
CCLE_transcripts_ADAR <- data.frame(CCLE_transcripts_ADAR, row.names = CCLE_transcripts_ADAR$transcript_id)

#remove transcrpt_ids and gene_id
CCLE_transcripts_ADAR$transcript_id <- NULL
CCLE_transcripts_ADAR$gene_id <- NULL

#transpose
CCLE_transcripts_ADAR <- as.data.frame(t(CCLE_transcripts_ADAR))

#split apart cell_line identifiers and add new columns for cell line and site
out <- strsplit(as.character(row.names(CCLE_transcripts_ADAR)), "_", fixed = TRUE)

out <- do.call(rbind, out)

out <- as.data.frame(out)

CCLE_transcripts_ADAR$cell_line <- tolower(as.character(out$V1))
CCLE_transcripts_ADAR$Site <- as.character(out$V2)

#keep only breast cell lines
CCLE_transcripts_ADAR <- subset(CCLE_transcripts_ADAR, Site == "BREAST")

#rename columns for p150 and p110 trnascripts
colnames(CCLE_transcripts_ADAR)[colnames(CCLE_transcripts_ADAR) == "ENST00000368474.4"] <- "p150"
colnames(CCLE_transcripts_ADAR)[colnames(CCLE_transcripts_ADAR) == "ENST00000368471.3"] <- "p110"


Dependency_CCLE_transcripts <- merge(CCLE_transcripts_ADAR, Dependency_CCLE, by = "cell_line")


#convert isoform expression to log2(TPM+1) to match gene level expression
ggplot(Dependency_CCLE_transcripts, aes(`DHX9 (ENSG00000135829)`, log2(p110+1))) + 
  geom_point() + theme_science() + 
  geom_smooth(method = "lm", colour = palette2[2]) + labs(y = "ADAR p110 Expression\n (log2(TPM+1))", x = "DHX9 Expression\n (log2(TPM+1))", colour = "") +
  theme(legend.position = "bottom") + stat_cor(method="pearson", alpha = 1)
ggsave("dhx9_p110.tiff", height = 3.7, width = 3.5)

ggplot(Dependency_CCLE_transcripts, aes(`DHX9 (ENSG00000135829)`, log2(p150+1))) + 
  geom_point() + theme_science() +
  geom_smooth(method = "lm", colour = palette2[2]) + labs(y = "ADAR p150 Expression\n (log2(TPM+1))", x = "DHX9 Expression\n (log2(TPM+1))", colour = "") +
  theme(legend.position = "bottom") + stat_cor(method="pearson", alpha = 1)
ggsave("dhx9_p150.tiff", height = 3.7, width = 3.5)




