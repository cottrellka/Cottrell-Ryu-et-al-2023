setwd("/Users/cottr/Box Sync/APEX-ADAR/")

library(ggplot2)
library(edgeR)
library(dplyr)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(data.table)
library(topGO)
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

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(DESeq2)
library(MetBrewer)
library(drawProteins)
library(gplots)
library(Vennerable)

library('biomaRt')
ensembl <- useMart('ensembl', dataset="hsapiens_gene_ensembl")


library(enrichR)

palette <- met.brewer(name="Veronese", n=7, type="discrete")



met.brewer(name="Veronese", n=7, type="discrete")

palette2 <- c(palette[7], palette[3])
palette3 <- c(palette[1], palette[6], palette[2])


SKBR3 <- read.delim("SKBR3_MS_Raw.txt")

SKBR3_short <- SKBR3

row.names(SKBR3_short) <- SKBR3_short$Accession.Number

SKBR3_short$Identified.Proteins..2630. <- NULL
SKBR3_short$Accession.Number <- NULL
SKBR3_short$Molecular.Weight <- NULL

SKBR3_p110 <- SKBR3_short

SKBR3_p110$APEX.p150_R1 <- NULL
SKBR3_p110$APEX.p150_R2 <- NULL

SKBR3_p110 <- as.matrix(SKBR3_p110)

head(SKBR3_p110)

coldata <- DataFrame(condition=factor(c("1","1","2","2")), 
                     row.names=as.character(colnames(SKBR3_p110)))
coldata

ddsHTSeq <- DESeqDataSetFromMatrix(countData= SKBR3_p110, 
                                   colData = coldata, 
                                   design = ~ condition )

ddsHTSeq

dds <- DESeq(ddsHTSeq)

res_skbr3 <- results(dds)

res_SKBR3_Norm <- lfcShrink(dds, coef=2, type="normal")

plotMA(res_skbr3)
plotMA(res_SKBR3_Norm)


RNA_deseq_summary_p110 <- as.data.frame(summary(res_skbr3))

sum(res_skbr3$padj <0.01, na.rm=TRUE)

res_skbr3_p110 <- as.data.frame(res_skbr3)
res_skbr3_p110_norm <- as.data.frame(res_SKBR3_Norm)

res_skbr3_p110$gene <- row.names(res_skbr3_p110)

res_skbr3_p110_norm$gene <- row.names(res_skbr3_p110_norm)

res_skbr3_p110$Sig <- ifelse(res_skbr3_p110$padj<0.05, "<0.05", ">0.05")
write.csv(res_skbr3_p110, "SKBR3_deseq_p110.csv")

res_skbr3_p110_norm$Sig <- ifelse(res_skbr3_p110_norm$padj<0.05, "<0.05", ">0.05")

ggplot(res_skbr3_p110, aes(log2FoldChange, -log10(padj), colour=Sig)) + geom_point(alpha=0.2) +
  theme_science() + scale_colour_manual(values=cbPalette) + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="Green") +
  annotate('rect', xmin = -0.95, xmax = 0.995, ymin = -0.1, ymax = Inf,alpha = .75, fill="white")
ggsave("volcano_SKBR3.tiff", height=4, width=4, unit="in")


ggplot(res_skbr3_p110_norm, aes(log2FoldChange, -log10(padj), colour=Sig)) + geom_point(alpha=0.2) +
  theme_science() + scale_colour_manual(values=cbPalette) + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="Green") +
  annotate('rect', xmin = -0.95, xmax = 0.995, ymin = -0.1, ymax = Inf,alpha = .75, fill="white")

res_skbr3_p110_up <- subset(res_skbr3_p110, Sig == "<0.05" & log2FoldChange > 0.5)
res_skbr3_p110_norm_up <- subset(res_skbr3_p110_norm, Sig == "<0.05" & log2FoldChange > 0.5)

write.csv(res_skbr3_p110_up, "deseq_SKBR3_up.csv")


SKBR3_p110 <- SKBR3_short

SKBR3_p110$APEX.p150_R1 <- NULL
SKBR3_p110$APEX.p150_R2 <- NULL

SKBR3_up_complete <- merge(res_skbr3_p110_up, SKBR3_p110, by = "row.names")

write.csv(SKBR3_up_complete, "deseq_p110_up_skbr3_complete.csv")


ddsHTSeq <- NULL
coldata <- NULL


MCF7 <- read.delim("MCF7_MS_Raw.txt")

MCF7_short <- MCF7

row.names(MCF7_short) <- MCF7_short$Accession.Number

MCF7_short$Identified.Proteins..3198. <- NULL
MCF7_short$Accession.Number <- NULL
MCF7_short$Molecular.Weight <- NULL

MCF7_p110 <- MCF7_short

MCF7_p110$APEX.p150_R1 <- NULL
MCF7_p110$APEX.p150_R2 <- NULL

MCF7_p110 <- as.matrix(MCF7_p110)

head(MCF7_p110)

coldata <- DataFrame(condition=factor(c("1","2","1","2")), 
                     row.names=as.character(colnames(MCF7_p110)))
coldata

ddsHTSeq <- DESeqDataSetFromMatrix(countData= MCF7_p110, 
                                   colData = coldata, 
                                   design = ~ condition )

ddsHTSeq

dds <- DESeq(ddsHTSeq)

res_MCF7 <- results(dds)

res_MCF7_Norm <- lfcShrink(dds, coef=2, type="normal")

plotMA(res_MCF7)
plotMA(res_MCF7_Norm)


RNA_deseq_summary_p110 <- as.data.frame(summary(res_MCF7))

sum(res_MCF7$padj <0.01, na.rm=TRUE)

res_MCF7_p110 <- as.data.frame(res_MCF7)
res_MCF7_p110_norm <- as.data.frame(res_MCF7_Norm)

res_MCF7_p110$gene <- row.names(res_MCF7_p110)

res_MCF7_p110_norm$gene <- row.names(res_MCF7_p110_norm)

res_MCF7_p110$Sig <- ifelse(res_MCF7_p110$padj<0.05, "<0.05", ">0.05")
write.csv(res_MCF7_p110, "MCF7_deseq_p110.csv")

res_MCF7_p110_norm$Sig <- ifelse(res_MCF7_p110_norm$padj<0.05, "<0.05", ">0.05")

ggplot(res_MCF7_p110, aes(log2FoldChange, -log10(padj), colour=Sig)) + geom_point(alpha=0.2) +
  theme_science() + scale_colour_manual(values=cbPalette) + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="Green") +
  annotate('rect', xmin = -0.95, xmax = 0.995, ymin = -0.1, ymax = Inf,alpha = .75, fill="white")
ggsave("volcano_MCF7.tiff", height=4, width=4, unit="in")


ggplot(res_MCF7_p110_norm, aes(log2FoldChange, -log10(padj), colour=Sig)) + geom_point(alpha=0.2) +
  theme_science() + scale_colour_manual(values=cbPalette) + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="Green") +
  annotate('rect', xmin = -0.95, xmax = 0.995, ymin = -0.1, ymax = Inf,alpha = .75, fill="white")

res_MCF7_p110_up <- subset(res_MCF7_p110, Sig == "<0.05" & log2FoldChange > 0.5)
res_MCF7_p110_norm_up <- subset(res_MCF7_p110_norm, Sig == "<0.05" & log2FoldChange > 0.5)

write.csv(res_MCF7_p110_up, "deseq_MCF7_up.csv")




ddsHTSeq <- NULL
coldata <- NULL


MB231 <- read.delim("MB231_MS_raw.txt")

MB231_short <- MB231

MB231$Accession.Number

row.names(MB231_short) <- MB231_short$Accession.Number

MB231_short$Identified.Proteins..2097. <- NULL
MB231_short$Accession.Number <- NULL
MB231_short$Molecular.Weight <- NULL

MB231_p110 <- MB231_short

MB231_p110 <- as.matrix(MB231_p110)

head(MB231_p110)

coldata <- DataFrame(condition=factor(c("1","1","2","2")),
                      row.names=as.character(colnames(MB231_p110)))
coldata

ddsHTSeq_mb <- DESeqDataSetFromMatrix(countData= MB231_p110,
                                    colData = coldata,
                                    design = ~ condition )

ddsHTSeq_mb

dds_mb <- DESeq(ddsHTSeq_mb)

res_mb <- results(dds_mb)

res_mb

RNA_deseq_summary_p110 <- as.data.frame(summary(res_mb))

sum(res_mb$padj <0.01, na.rm=TRUE)

res_p110_mb <- as.data.frame(res_mb)

res_p110_mb$gene <- row.names(res_p110_mb)

res_p110_mb$Sig <- ifelse(res_p110_mb$padj<0.05, "<0.05", ">0.05")
write.csv(res_p110_mb, "deseq_p110_MB231.csv")

ggplot(res_p110_mb, aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha = 0.2) +
   theme_science() + scale_colour_manual(values=cbPalette) + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") +
   geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="Green") +
   annotate('rect', xmin = -0.995, xmax = 0.995, ymin = -0.5, ymax = Inf,alpha = .75, fill="white")
ggsave("volcano_p110_MB231.tiff", height=4, width=4, unit="in")

res_p110_mb_up <- subset(res_p110_mb, Sig == "<0.05" & log2FoldChange > 0.5)

write.csv(res_p110_mb_up, "deseq_p110_up_MB231.csv")


res_p110_mb$cell_line <- "MDA-MB-231"
res_MCF7_p110$cell_line <- "MCF-7"
res_skbr3_p110$cell_line <- "SK-BR-3"

res_all <- rbind(res_p110_mb, res_MCF7_p110, res_skbr3_p110)

ggplot(res_all, aes(log2FoldChange, -log10(padj), colour = cell_line)) + geom_point(alpha = 0.5) +
  theme_science() + scale_colour_manual(values=palette3) + labs(x = "Log2 Fold Change", y = "-log10 (FDR)", colour = "") + 
  annotate('rect', xmin = -1, xmax = 1, ymin = -0.5, ymax = Inf,alpha = .75, fill="white") +
  annotate('rect', xmin = -10.2, xmax = 10, ymin = -0.5, ymax = -log10(0.05),alpha = .75, fill="white") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="black") +
  geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") +
  theme(legend.position = "bottom")
ggsave("volcano_p110.tiff", height=4, width=4, unit="in")


write.csv(res_all, "apex_results.csv")


out <- strsplit(res_MCF7_p110_up$gene, "|", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

out_2 <- strsplit(as.character(out$V3), "_", fixed = TRUE)
out_2 <- do.call(rbind, out_2)
out_2 <- as.data.frame(out_2)

res_MCF7_p110_up$symbol <- out_2$V1

res_MCF7_p110_up$protein <- out$V2


out <- strsplit(res_skbr3_p110_up$gene, "|", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

out_2 <- strsplit(as.character(out$V3), "_", fixed = TRUE)
out_2 <- do.call(rbind, out_2)
out_2 <- as.data.frame(out_2)

res_skbr3_p110_up$symbol <- out_2$V1

res_skbr3_p110_up$protein <- out$V2


out <- strsplit(res_p110_mb_up$gene, "|", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

out_2 <- strsplit(as.character(out$V3), "_", fixed = TRUE)
out_2 <- do.call(rbind, out_2)
out_2 <- as.data.frame(out_2)

res_p110_mb_up$symbol <- out_2$V1

res_p110_mb_up$protein <- out$V2

list_venn <- list(`SKBR3` = res_skbr3_p110_up$symbol,
                  `MCF-7` = res_MCF7_p110_up$symbol,
                  `MDA-MB-231` = res_p110_mb_up$symbol)

Vstem <- venn(list_venn)

plot(Vstem, doWeights = TRUE, show = list(Faces = TRUE, DarkMatter = FALSE), doEuler=TRUE)

test <- Reduce(intersect, list_venn)

w <- compute.Venn(Venn(Sets=list_venn))
gp <- VennThemes(w)
gp[["Set"]][["Set1"]]$col <-  "#99610a"
gp[["Set"]][["Set2"]]$col <-  "#67322e"
gp[["Set"]][["Set3"]]$col <-  "#175449"
gp[["Face"]][["001"]]$fill <- "white"
gp[["Face"]][["010"]]$fill <- "white"
gp[["Face"]][["011"]]$fill <- "grey"
gp[["Face"]][["111"]]$fill <- "grey"
gp[["Face"]][["101"]]$fill <- "grey"
gp[["Face"]][["110"]]$fill <- "grey"
gp[["Face"]][["100"]]$fill <- "white"
gp[["SetText"]][["Set1"]]$col <- "#99610a"
gp[["SetText"]][["Set2"]]$col <- "#67322e"
gp[["SetText"]][["Set3"]]$col <- "#175449"
grid.newpage()
tiff('venn_cell_lines.tiff', res = 600, width = 3500, height = 3500)
plot(w, gp = gp)
dev.off()


all_up <- rbind(res_skbr3_p110_up, res_MCF7_p110_up, res_p110_mb_up)

write.csv(all_up, "all_up.csv")



subcell <- fread("SubCellBarcode.MCF7.csv")

subcell_nuclear <- subset(subcell, Neighborhood.Class == "Nuclear")

list_venn <- list(`APEX2-ADAR` = all_up$symbol,
                  Nuclear = subcell_nuclear$Protein)

Vstem <- Venn(list_venn)

plot(Vstem, doWeights = TRUE, show = list(Faces = TRUE, DarkMatter = FALSE), doEuler=TRUE)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

test <- Reduce(intersect, list_venn)

w <- compute.Venn(Venn(Sets=list_venn))
gp <- VennThemes(w)
gp[["Set"]][["Set1"]]$col <-  "Black"
gp[["Set"]][["Set2"]]$col <-  palette[3]
gp[["Face"]][["01"]]$fill <- "white"
gp[["Face"]][["10"]]$fill <- "white"
gp[["Face"]][["11"]]$fill <- "Grey"
gp[["SetText"]][["Set1"]]$col <- "Black"
gp[["SetText"]][["Set2"]]$col <- palette[3]
grid.newpage()
tiff('venn_nuclear.tiff', res = 600, width = 3700, height = 3500)
plot(w, gp = gp)
dev.off()

subcell_n3 <- subset(subcell, Compartment.Class == "N3")

list_venn <- list(`APEX2-ADAR` = all_up$symbol,
                  Nucleolar = subcell_n3$Protein)

Vstem <- Venn(list_venn)

plot(Vstem, doWeights = TRUE, show = list(Faces = TRUE, DarkMatter = FALSE), doEuler=TRUE)

test <- Reduce(intersect, list_venn)

w <- compute.Venn(Venn(Sets=list_venn))
gp <- VennThemes(w)
gp[["Set"]][["Set1"]]$col <-  "Black"
gp[["Set"]][["Set2"]]$col <-  palette[6]
gp[["Face"]][["01"]]$fill <- "white"
gp[["Face"]][["10"]]$fill <- "white"
gp[["Face"]][["11"]]$fill <- "Grey"
gp[["SetText"]][["Set1"]]$col <- "Black"
gp[["SetText"]][["Set2"]]$col <- palette[6]
grid.newpage()
tiff('venn_nucleolar.tiff', res = 600, width = 3700, height = 3500)
plot(w, gp = gp)
dev.off()





dbs <- listEnrichrDbs()

genes <- unique(all_up$protein)

annotation <- getBM(attributes=c("uniprotswissprot", "hgnc_symbol"), 
                    filters="uniprotswissprot", values= genes, mart=ensembl)

all_up <- merge(all_up, annotation, by.x = "protein", by.y = "uniprotswissprot")

genes <- unique(all_up$hgnc_symbol)

enriched_biological <- enrichr(genes, "GO_Biological_Process_2021")

go_biological <- as.data.frame(enriched_biological)

go_biological <- subset(go_biological, GO_Biological_Process_2021.Adjusted.P.value < 0.05)

write.csv(go_biological, "go_biological.csv")

out <- strsplit(go_biological$GO_Biological_Process_2021.Overlap, "/", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

go_biological$number <- out$V1


go_biological_top <- dplyr::arrange(go_biological, desc(-GO_Biological_Process_2021.Adjusted.P.value))

go_biological_top <- go_biological_top[1:30,]

go_biological_top$GO_Biological_Process_2021.Term <- gsub("\\(.*", "", go_biological_top$GO_Biological_Process_2021.Term)


ggplot(go_biological_top,
       aes(-log10(GO_Biological_Process_2021.Adjusted.P.value), 
           reorder(GO_Biological_Process_2021.Term, -GO_Biological_Process_2021.Adjusted.P.value), 
           size = as.numeric(number), colour = as.numeric(GO_Biological_Process_2021.Odds.Ratio))) + 
  geom_point() + theme_science() + labs(x = "-log10(FDR)", y = "GO Biological Process", colour = "Odds Ratio", size = "Proteins \nIdentified") +
  scale_colour_continuous(type = "viridis")
ggsave("go_biological.tiff", height = 6, width = 15, units = "in")


enriched_molecular <- enrichr(genes, "GO_Molecular_Function_2021")

go_molecular <- as.data.frame(enriched_molecular)

go_molecular <- subset(go_molecular, GO_Molecular_Function_2021.Adjusted.P.value < 0.05)

write.csv(go_molecular, "go_molecular.csv")

out <- strsplit(go_molecular$GO_Molecular_Function_2021.Overlap, "/", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

go_molecular$number <- out$V1


go_molecular_top <- dplyr::arrange(go_molecular, desc(-GO_Molecular_Function_2021.Adjusted.P.value))

go_molecular_top <- go_molecular_top[1:22,]

go_molecular_top$GO_Molecular_Function_2021.Term <- gsub("\\(.*", "", go_molecular_top$GO_Molecular_Function_2021.Term)


ggplot(go_molecular_top,
       aes(-log10(GO_Molecular_Function_2021.Adjusted.P.value), 
           reorder(GO_Molecular_Function_2021.Term, -GO_Molecular_Function_2021.Adjusted.P.value), 
           size = as.numeric(number), colour = as.numeric(GO_Molecular_Function_2021.Odds.Ratio))) + 
  geom_point() + theme_science() + labs(x = "-log10(FDR)", y = "GO molecular Process", colour = "Odds Ratio", size = "Proteins \nIdentified") +
  scale_colour_continuous(type = "viridis")
ggsave("go_molecular.tiff", height = 6, width = 15, units = "in")




enriched_cellular <- enrichr(genes, "GO_Cellular_Component_2021")

go_cellular <- as.data.frame(enriched_cellular)

go_cellular <- subset(go_cellular, GO_Cellular_Component_2021.Adjusted.P.value < 0.05)

write.csv(go_cellular, "go_cellular.csv")

out <- strsplit(go_cellular$GO_Cellular_Component_2021.Overlap, "/", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

go_cellular$number <- out$V1


go_cellular_top <- dplyr::arrange(go_cellular, desc(-GO_Cellular_Component_2021.Adjusted.P.value))

go_cellular_top$GO_Cellular_Component_2021.Term <- gsub("\\(.*", "", go_cellular_top$GO_Cellular_Component_2021.Term)


ggplot(go_cellular_top,
       aes(-log10(GO_Cellular_Component_2021.Adjusted.P.value), 
           reorder(GO_Cellular_Component_2021.Term, -GO_Cellular_Component_2021.Adjusted.P.value), 
           size = as.numeric(number), colour = as.numeric(GO_Cellular_Component_2021.Odds.Ratio))) + 
  geom_point() + theme_science() + labs(x = "-log10(FDR)", y = "GO cellular Process", colour = "Odds Ratio", size = "Proteins \nIdentified") +
  scale_colour_continuous(type = "viridis")
ggsave("go_cellular.tiff", height = 6, width = 15, units = "in")



go_biological_top10 <- go_biological_top[1:10,]
go_biological_top10$type <- 1
colnames(go_biological_top10)[1] <- "GO_term"
colnames(go_biological_top10) <- gsub("GO_Biological_Process.", "", colnames(go_biological_top10))

go_molecular_top10 <- go_molecular_top[1:10,]
go_molecular_top10$type <- 2
colnames(go_molecular_top10)[1] <- "GO_term"
colnames(go_molecular_top10) <- gsub("GO_Molecular_Function.", "", colnames(go_molecular_top10))

go_cellular_top10 <- go_cellular_top[1:10,]
go_cellular_top10$type <- 3
colnames(go_cellular_top10)[1] <- "GO_term"
colnames(go_cellular_top10) <- gsub("GO_Cellular_Component.", "", colnames(go_cellular_top10))

go_top <- rbind(go_biological_top10, go_molecular_top10, go_cellular_top10)

go_top$GO_term <- as.factor(go_top$GO_term)

go_top$GO_term <- factor(go_top$GO_term, levels=go_top$GO_term[order(-go_top$`2021.Adjusted.P.value`)], ordered = TRUE)

ggplot(go_top, aes(-log10(`2021.Adjusted.P.value`), reorder(`GO_term`, - as.numeric(type)), size = as.numeric(`number`), colour = as.numeric(`2021.Odds.Ratio`))) + 
  geom_point() + theme_science() + labs(x = "-log10(FDR)", y = "GO Cellular Component    GO Molecular Function   GO Biological Process", colour = "Fold \nEnrichment", size = "Proteins \nIdentified") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + geom_hline(yintercept = c(10.5,20.5,30.5)) + 
  theme(legend.position = c(0.3,-0.05), legend.direction = "horizontal", plot.margin = grid::unit(c(0.1,0.1,5,0.1), 'lines'))
ggsave("go_all.tiff", height = 8, width = 7.5, units = "in")

ggplot(go_top, aes(-log10(`2021.Adjusted.P.value`), reorder(`GO_term`, - as.numeric(type)), size = as.numeric(`number`), colour = as.numeric(`2021.Odds.Ratio`))) + 
  geom_point() + theme_science() + labs(x = "-log10(FDR)", y = "GO Cellular Component    GO Molecular Function   GO Biological Process", colour = "Fold \nEnrichment", size = "Proteins \nIdentified") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + geom_hline(yintercept = c(10.5,20.5,30.5)) + 
  theme(legend.position = "bottom", legend.justification = c(1,0), axis.title.y = element_text(vjust = 0))
ggsave("go_all.tiff", height = 8, width = 7.5, units = "in")



all_up <- read.csv("all_up.csv")

all_helicases <- all_up[grepl("DDX|DHX", all_up$symbol),]

helicases <- unique(all_helicases$protein)


prot_data <- drawProteins::get_features( "Q08211 Q92841 P17844 Q96GQ7 Q9NR30 Q9GZR7 Q9NVP1 Q8TDD1 Q9BQ39 P19525 P55265")

drawProteins::feature_to_dataframe(prot_data) -> prot_data


p <- draw_canvas(prot_data)
p <- draw_chains(p, prot_data)
p <- draw_domains(p, prot_data)

new_order <- function(x) { 
  if(grepl("E2AK2", x)) y <- 11
  if(grepl("DSRAD", x)) y <- 10
  if(grepl("DHX9", x)) y <- 9
  if(grepl("DDX5", x)) y <- 8
  if(grepl("DDX17", x)) y <- 7
  if(grepl("DDX18", x)) y <- 6
  if(grepl("DDX21", x)) y <- 5
  if(grepl("DDX24", x)) y <- 4 
  if(grepl("DDX27", x)) y <- 3
  if(grepl("DDX50", x)) y <- 2
  if(grepl("DDX54", x)) y <- 1  
  return(y)
}

#apply stop codon function
prot_data$order <- sapply(prot_data$entryName, new_order)

prot_data$entryName <- gsub("_HUMAN", "", prot_data$entryName)

prot_data$entryName <- gsub("DSRAD", "ADAR1", prot_data$entryName)

prot_data$entryName <- gsub("E2AK2", "PKR", prot_data$entryName)

prot_data$description <- gsub("M 1|M 2|M 3", "M", prot_data$description)

prot_data$description <- gsub("g 1|g 2", "g", prot_data$description)

prot_data$description <- gsub("DRBM", "dsRBD", prot_data$description)

p <- draw_canvas(prot_data)
p <- draw_chains(p, prot_data, fill = "lightgrey", outline = "white")
p <- draw_domains(p, prot_data, label_domains = FALSE) 

p + theme_science() + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), 
                            axis.ticks = element_blank(),
                            axis.text = element_blank(), axis.title.x = element_blank(), legend.position = "top") +
  labs(fill = "") + scale_fill_manual(values = c(palette[2], palette[1], palette[4], palette[5], palette[7], palette[3]))
ggsave("domains.tiff", width = 10, height = 4, units  = "in")
ggsave("domains_tall.tiff", width = 10, height = 10, units  = "in")

# prot_data <- drawProteins::get_features( "Q08211 Q92841 P17844 Q96GQ7 Q9NR30 Q9GZR7 Q9NVP1 Q8TDD1 Q9BQ39 P19525 P55265")
# 
# drawProteins::feature_to_dataframe(prot_data) -> prot_data
# 
# new_order <- function(x) { 
#   if(grepl("E2AK2", x)) y <- 22
#   if(grepl("DSRAD", x)) y <- 20
#   if(grepl("DHX9", x)) y <- 18
#   if(grepl("DDX5", x)) y <- 16
#   if(grepl("DDX17", x)) y <- 14
#   if(grepl("DDX18", x)) y <- 12
#   if(grepl("DDX21", x)) y <- 10
#   if(grepl("DDX24", x)) y <- 8 
#   if(grepl("DDX27", x)) y <- 6
#   if(grepl("DDX50", x)) y <- 4
#   if(grepl("DDX54", x)) y <- 2  
#   return(y)
# }
# 
# #apply stop codon function
# prot_data$order <- sapply(prot_data$entryName, new_order)
# 
# prot_data$entryName <- gsub("_HUMAN", "", prot_data$entryName)
# 
# prot_data$entryName <- gsub("DSRAD", "ADAR1", prot_data$entryName)
# 
# prot_data$entryName <- gsub("E2AK2", "PKR", prot_data$entryName)
# 
# prot_data$description <- gsub("M 1|M 2|M 3", "M", prot_data$description)
# 
# prot_data$description <- gsub("g 1|g 2", "g", prot_data$description)
# 
# prot_data$description <- gsub("DRBM", "dsRBD", prot_data$description)
# 
# p <- draw_canvas(prot_data)
# p <- draw_chains(p, prot_data, fill = "lightgrey", outline = "white")
# p <- draw_domains(p, prot_data, label_domains = FALSE) 
# 
# p + theme_science() + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), 
#                             axis.ticks = element_blank(),
#                             axis.text = element_blank(), axis.title.x = element_blank(), legend.position = "top") +
#   labs(fill = "") + scale_fill_manual(values = c(palette[2], palette[1], palette[4], palette[5], palette[7], palette[3]))
# ggsave("domains.tiff", width = 10, height = 4, units  = "in")
# ggsave("domains_tall.tiff", width = 10, height = 10, units  = "in")






prot_data_dhx9_adar <- subset(prot_data, entryName %in% c("ADAR", "DHX9"))

prot_data_dhx9_adar <- subset(prot_data_dhx9_adar, type %in% c("CHAIN", "DOMAIN"))

#write.csv(prot_data_dhx9_adar, "prot_data_dhx9_adar.csv")

prot_data_dhx9_adar <- read.csv("prot_data_dhx9_adar.csv")

prot_data_dhx9_adar$description <- gsub("DRBM", "dsRBD", prot_data_dhx9_adar$description)



p2 <- draw_canvas(prot_data_dhx9_adar)
p2 <- draw_chains(p2, prot_data_dhx9_adar, fill = "lightgrey", outline = "white")
p2 <- draw_domains(p2, prot_data_dhx9_adar, label_domains = FALSE) 

p2 + theme_science() + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), 
                            axis.ticks = element_blank(), axis.text = element_blank() ,
                            axis.title.x = element_blank(), legend.position = "top") +
  labs(fill = "") + scale_fill_manual(values = c(palette[2], palette[1], palette[7], palette[4], palette[5], palette[3], "black")) +
  annotate("text", x = 417, y = 3, label = "K417R") + 
  annotate("segment", x = 417, xend = 417, y = 1.75, yend = 2.7, colour = "black")
ggsave("domains_2.tiff", width = 6, height = 4, units  = "in")
ggsave("domains_2_wide.tiff", width = 10, height = 4, units  = "in")



