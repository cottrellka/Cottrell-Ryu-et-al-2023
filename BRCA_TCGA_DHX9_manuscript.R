library(RTCGA)
library(Biobase)
library(affy)
library(genefilter)
library(XML)
library(reshape2)
library(survival)
library(survminer)
require("survival")
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(TCGAretriever)
library(extrafont)
library(DESeq2)
library(tidyverse)
library(edgeR)
library(DescTools)
library(MetBrewer)
library(ggpubr)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

palette <- met.brewer(name="Veronese", n=7, type="discrete")

met.brewer(name="Veronese", n=7, type="discrete")

palette2 <- c(palette[7], palette[3])

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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

setwd("/Users/cottr/Box Sync/TCGA/")


# read in RNAseq data for breast cancer samples from TCGA - 
RNAseq_BRCA <- fread("/Users/cottr/Box Sync/TCGA/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt")

#first row has relevant column names, make col.names match the first row
colnames(RNAseq_BRCA) <- paste(colnames(RNAseq_BRCA), RNAseq_BRCA[1,], sep= "-")

#hybridization reference has gene name and id as a string, split it and make a data.frame of the split
out <- strsplit(as.character(RNAseq_BRCA$`Hybridization REF`), "\\|")
test <- as.data.frame(do.call(rbind, out))


#get just the columns with raw counts
RNAseq_BRCA <- dplyr::select(RNAseq_BRCA, contains("raw_count"))

#Make a new column for gene name using the split from above
RNAseq_BRCA$Gene  <- test$V1

#remove the first row
RNAseq_BRCA <- RNAseq_BRCA[-1,]

#get rid of genes with ? or duplicate names
RNAseq_BRCA <- subset(RNAseq_BRCA, Gene != "?")
RNAseq_BRCA <- subset(RNAseq_BRCA, Gene != "SLC35E2")

#make a new data.frame with RNAseq_BRCA data with row.names = gene name
RNAseq_BRCA <- data.frame(RNAseq_BRCA, row.names = RNAseq_BRCA$Gene)

#get rid of the the gene name column
RNAseq_BRCA$Gene <- NULL

#make a new data.frame for RNAseq_BRCA where all of the raw counts are numeric
RNAseq_BRCA_df <- as.data.frame(row.names = row.names(RNAseq_BRCA), lapply(RNAseq_BRCA, function(x) as.numeric(as.character(x))))


str(RNAseq_BRCA_df)

row.names(RNAseq_BRCA_df)

#the next few lines are based off of this example: https://www.biostars.org/p/153013/
#quick check to see how many normal and primary samples there are. 
#In TCGA barcode Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29.
#the 14th character will be 0 for tumor and 1 for normal
table(substr(colnames(RNAseq_BRCA_df),14,14))

#index tumor (t) or normal (n) based on barcode value
n_index <- which(substr(colnames(RNAseq_BRCA_df),14,14) == '1')
t_index <- which(substr(colnames(RNAseq_BRCA_df),14,14) == '0')

#get CPM values for RNAseq data
logCPM <- cpm(RNAseq_BRCA_df, log=TRUE)

#define scal function, calculates a modified z-score, see link above for more details
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}

makeStars <- function(x){
  stars <- c("***", "**", "*", "ns")
  vec <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

#apply scal function to the logCPM data 
z_rna <- scal(logCPM,logCPM[,n_index])

#transpose z_rna and make new RNAseq_BRCA data.frame
RNAseq_BRCA <- as.data.frame(t(z_rna))


#make a column for sample barcode based on row.names
RNAseq_BRCA$Sample <- row.names(RNAseq_BRCA)

#split the sample barcode
out <- strsplit(as.character(RNAseq_BRCA$Sample), "\\.")
test <- as.data.frame(do.call(rbind, out))


#make the patient barcode using the first three values of the split
test$bcr_patient_barcode  <- do.call(paste, c(test[c("V1", "V2", "V3")], sep = "-")) 

#add the patient_barcode to the data.frame
RNAseq_BRCA$patient_barcode <- test$bcr_patient_barcode

#make the sample barcode using the first three values of the split
test$Sample_Barcode <- do.call(paste, c(test[c("V1", "V2", "V3","V4")], sep = "-")) 

#add the sample_barcode to the data.frame
RNAseq_BRCA$Sample_Barcode <- test$Sample_Barcode

#add the sample type
RNAseq_BRCA$Sample_Type <- test$V4

#rename the sample types, add DROP for the B samples for each normal, primary or metastatic
RNAseq_BRCA$Sample_Type <- gsub("11A", "Normal", RNAseq_BRCA$Sample_Type)
RNAseq_BRCA$Sample_Type <- gsub("11B", "DROP", RNAseq_BRCA$Sample_Type)

RNAseq_BRCA$Sample_Type <- gsub("01A", "Primary", RNAseq_BRCA$Sample_Type)
RNAseq_BRCA$Sample_Type <- gsub("01B", "DROP", RNAseq_BRCA$Sample_Type)

RNAseq_BRCA$Sample_Type <- gsub("06A", "Metastatic", RNAseq_BRCA$Sample_Type)
RNAseq_BRCA$Sample_Type <- gsub("06B", "DROP", RNAseq_BRCA$Sample_Type)

#drop the b samples
RNAseq_BRCA <- subset(RNAseq_BRCA, !Sample_Type == "DROP")

#subset only the primary tumor samples
RNAseq_BRCA_primary <- subset(RNAseq_BRCA, Sample_Type == "Primary")

RNAseq_BRCA_primary$Sample_Barcode

#load annotation file from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4911051/ 
Annotation <- fread("journal.pone.0157368.s008.txt")

#keep only columns 2-4 of annotation file
Annotation <- Annotation[,2:4]

#merge with RNAseq data, patient barcode works here because there are no normal or metastic samples
RNAseq_BRCA_primary <- merge(RNAseq_BRCA_primary, Annotation, by.y = "BARCODE", by.x = "patient_barcode", all.x = TRUE)

#subset only the normal samples
RNAseq_BRCA_normal <- subset(RNAseq_BRCA, Sample_Type == "Normal")

#defne TNBC and PAM50 as Normal or Normal Breast
RNAseq_BRCA_normal$TNBC <- "Normal"
RNAseq_BRCA_normal$PAM50 <- "Normal Breast"

#bind primary and normal data
RNAseq_BRCA_primary_normal <- rbind(RNAseq_BRCA_normal, RNAseq_BRCA_primary)

#make plots of interest and perform hypothesis testing
aov_DHX9 <- aov(DHX9 ~ TNBC, RNAseq_BRCA_primary_normal)

summary(aov_DHX9)

RNAseq_BRCA_primary_normal$TNBC <- factor(as.factor(RNAseq_BRCA_primary_normal$TNBC), levels = c("Normal", "NO", "YES"))

DT_tnbc <- DunnettTest(RNAseq_BRCA_primary_normal$DHX9, RNAseq_BRCA_primary_normal$TNBC)

DT_tnbc <- as.data.frame(DT_tnbc$Normal[,1:4])

out <- strsplit(as.character(row.names(DT_tnbc)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_tnbc$group1 <- out$V1
DT_tnbc$group2 <- out$V2
DT_tnbc$labels <- makeStars(DT_tnbc$pval)
DT_tnbc$y.position <- max(RNAseq_BRCA_primary_normal$DHX9*1.1)


ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(TNBC)), aes(as.factor(TNBC), DHX9)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "NO", "YES"), labels=c("NO" = "Non TNBC", "YES" = "TNBC")) + ylab("DHX9 Expression\n Z-Score") +
  geom_bracket(xmin = DT_tnbc$group2, 
               xmax = DT_tnbc$group1, 
               y.position = DT_tnbc$y.position, label = DT_tnbc$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16) 
ggsave("TNBC_v_DHX9.tiff", width = 5, height = 4, units = "in")

#make plots of interest and perform hypothesis testing
aov_DHX9 <- aov(DHX9 ~ PAM50, RNAseq_BRCA_primary_normal)

summary(aov_DHX9)

RNAseq_BRCA_primary_normal$PAM50 <- factor(as.factor(RNAseq_BRCA_primary_normal$PAM50), 
                                           levels = c("Normal Breast", "Normal", "Basal", "Her2", "LumA", "LumB"))

DTpam50 <- DunnettTest(RNAseq_BRCA_primary_normal$DHX9, RNAseq_BRCA_primary_normal$PAM50)



DT_pam50 <- DunnettTest(RNAseq_BRCA_primary_normal$DHX9, RNAseq_BRCA_primary_normal$PAM50)

DT_pam50 <- as.data.frame(DT_pam50$Normal[,1:4])

out <- strsplit(as.character(row.names(DT_pam50)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_pam50$group1 <- out$V1
DT_pam50$group2 <- out$V2
DT_pam50$labels <- makeStars(DT_pam50$pval)
DT_pam50$y.position <- max(RNAseq_BRCA_primary_normal$DHX9*1.1)


ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(PAM50)), aes(as.factor(PAM50), DHX9)) + xlab("") + 
  scale_x_discrete(limits= c("Normal Breast", "Normal", "Basal", "Her2", "LumA", "LumB"), labels=c("Normal\nBreast", "PAM50\nNormal", "Basal", "Her2", "LumA", "LumB")) +
  ylab("DHX9 Expression\n Z-Score") +
  geom_bracket(xmin = DT_pam50$group2, 
               xmax = DT_pam50$group1, 
               y.position = DT_pam50$y.position, label = DT_pam50$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16) 
ggsave("pam50_v_DHX9.tiff", width = 6, height = 4, units = "in")


#make plots of interest and perform hypothesis testing
aov_DHX9 <- aov(DHX9 ~ Sample_Type, RNAseq_BRCA_primary_normal)

summary(aov_DHX9)

RNAseq_BRCA$Sample_Type <- factor(as.factor(RNAseq_BRCA$Sample_Type), levels = c("Normal", "Primary", "Metastatic"))

DT_sample <- DunnettTest(RNAseq_BRCA$DHX9, RNAseq_BRCA$Sample_Type)

DT_sample <- as.data.frame(DT_sample$Normal)

out <- strsplit(as.character(row.names(DT_sample)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_sample$group1 <- out$V1
DT_sample$group2 <- out$V2
DT_sample$labels <- makeStars(DT_sample$pval)
DT_sample$y.position <- max(RNAseq_BRCA_primary_normal$DHX9*1.1)


ggplot(subset(RNAseq_BRCA, !is.na(Sample_Type)), aes(as.factor(Sample_Type), DHX9)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "Primary", "Metastatic")) + ylab("DHX9 Expression\n Z-Score") +
  geom_bracket(xmin = DT_sample$group2, 
               xmax = DT_sample$group1, 
               y.position = DT_sample$y.position, label = DT_sample$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) + geom_jitter(width = 0.2, alpha = 0.2) + theme_science(base_size = 16) 
ggsave("sample_v_DHX9.tiff", width = 5, height = 4, units = "in")




ggplot(RNAseq_BRCA_primary_normal, aes(DHX9, ADAR, colour = Sample_Type, alpha = Sample_Type)) + 
  geom_point() + theme_science() + scale_alpha_manual(values = c(0.5, 0.2), guide = "none") +
  geom_smooth(method = "lm") + labs(y = "ADAR Expression\n Z-score", x = "DHX9 Expression\n (Z-score)", colour = "") +
  scale_color_manual(values = palette2, labels = c("Normal Breast", "Primary BRCA")) + 
  theme(legend.position = "bottom") + stat_cor(method="pearson", alpha = 1)
ggsave("DHX9_ADAR.tiff", height = 4.2, width = 4)
cor.test(RNAseq_BRCA_isoforms$DHX9, RNAseq_BRCA_isoforms$p150)


#remove variables that aren't RNA expression
RNAseq_BRCA_primary <- RNAseq_BRCA_primary[,-20502:-20506]

#find correlaiton between ADAR expression and all other genes

CORS_rna_p <- matrix(nrow=20501, ncol=3) 
for (i in 2:length(RNAseq_BRCA_primary)) {
  a <- cor.test(RNAseq_BRCA_primary$ADAR, RNAseq_BRCA_primary[,i], method = "pearson")
  CORS_rna_p[i,] <- c(colnames(RNAseq_BRCA_primary)[i], a$estimate, a$p.value)
}

CORS_rna_p <- as.data.frame(na.omit(CORS_rna_p))

CORS_rna_p$CC <- as.numeric(as.character(CORS_rna_p$V2))
CORS_rna_p$p_value <- as.numeric(as.character(CORS_rna_p$V3))
CORS_rna_p$FDR <- p.adjust(CORS_rna_p$p_value, method = "fdr", n = length(CORS_rna_p$p_value))
CORS_rna_p$method = "Pearson"

#same as above but for spearman correlation
CORS_rna_s <- matrix(nrow=20501, ncol=3) 
for (i in 2:length(RNAseq_BRCA_primary)) {
  a <- cor.test(RNAseq_BRCA_primary$ADAR, RNAseq_BRCA_primary[,i], method = "spearman")
  CORS_rna_s[i,] <- c(colnames(RNAseq_BRCA_primary)[i], a$estimate, a$p.value)
}

CORS_rna_s <- as.data.frame(na.omit(CORS_rna_s))

CORS_rna_s$CC <- as.numeric(CORS_rna_s$V2)
CORS_rna_s$p_value <- as.numeric(CORS_rna_s$V3)

CORS_rna_s[CORS_rna_s == 0] <- 1E-20

CORS_rna_s$FDR <- p.adjust(CORS_rna_s$p_value, method = "fdr", n = length(CORS_rna_s$p_value))
CORS_rna_s$method = "Spearman"

CORS_rna_merge <- rbind(CORS_rna_p, CORS_rna_s)

#subset the correlation data for only DDX## and DHX## genes

CORS_rna_helicase <- subset(CORS_rna_merge, grepl("DDX|DHX", CORS_rna_merge$V1))

CORS_rna_helicase$rank <- CORS_rna_helicase$Pearson + CORS_rna_helicase$Spearman /2

out <- strsplit(CORS_rna_helicase$V1, " ", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

CORS_rna_helicase$Gene <- out$V1

#subset to only helicases identified by apex labeling, plot correlations

CORS_rna_helicase <- subset(CORS_rna_helicase, Gene %in% c("DHX9","DDX17", "DDX24", "DDX18", "DDX21", "DDX5", "DDX54", "DDX20", "DDX27"))

ggplot(CORS_rna_helicase, aes(CC, reorder(Gene, CC), size = -log10(FDR), colour = method)) + 
  geom_point() + theme_science() + labs(x = "Correlation Coefficient", y = "", colour = "") +
  scale_colour_manual(values = palette2)
ggsave("cors_tcga.tiff", height = 3, width = 4.5)

write.csv(CORS_rna_helicase, "CORS_rna_helicase.csv")






#read breast cancer clinical data
Clinical_BRCA <- readTCGA("/Users/cottr/Box Sync/TCGA/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.merged_only_clinical_clin_format.txt", dataType = 'clinical')

Clinical_BRCA$patient_barcode <- toupper(Clinical_BRCA$patient.bcr_patient_barcode)


#use the survivalTCGA function to extract survival data
BRCA_survival <- survivalTCGA(Clinical_BRCA)

BRCA_survival$times <- BRCA_survival$times/365

#merge survival data with primary tumor sample RNAseq
BRCA <- merge(BRCA_survival, RNAseq_BRCA_primary, by.x = "bcr_patient_barcode", by.y="patient_barcode")

#remove data with times <0
BRCA <- subset(BRCA, !times < 0)

#determine DHX9 expression cutoff using surv_cutpoint function

BRCA.surv_rnaseq.cut <- surv_cutpoint(
  BRCA,
  time = "times",
  event = "patient.vital_status",
  variables = "DHX9"
)
summary(BRCA.surv_rnaseq.cut)
plot(BRCA.surv_rnaseq.cut, "DHX9")

#score patients by DHX9 expression above or below cutoff
BRCA$DHX9_score <- ifelse(BRCA$DHX9 < BRCA.surv_rnaseq.cut$cutpoint$cutpoint, "low", "high")

#fit the survival curve and plot
fit_DHX9 <- survfit(Surv(times, patient.vital_status)
                    ~ DHX9_score , data = BRCA)


ggsurvplot(fit_DHX9, data = BRCA, size = 1,                 # change line size
           palette = palette2,# custom color palettes
           conf.int = TRUE,          # Add confidence interval
           pval = TRUE,              # Add p-value
           legend.title = "DHX9\nExpression",
           legend = "right",
           legend.labs = c("High", "Low"),
           xlab = "Years",
           risk.table = FALSE,        # Add risk table
           font.family = "Arial", 
           ggtheme = theme_science(base_size = 12), # Change ggplot2 theme
)
ggsave("DHX9_cutoff.tiff", width = 4, height = 3, units = "in")






#read isoform RNAseq data from FireBrowse
RNAseq_BRCA_isoforms <- fread("BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt")

#select only isoforms of interest (p150, p110, DHX9)
RNAseq_BRCA_isoforms <- subset(RNAseq_BRCA_isoforms, grepl("uc001ffh|uc001ffj|uc001ffi|uc001ffk|uc009wyd|uc001gpr", RNAseq_BRCA_isoforms$`Hybridization REF`))


#rename to protein name
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc001ffh.2", "p150-1", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc001ffi.2", "p150-2", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc001ffj.2", "p150-3", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc001ffk.2", "p110", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc009wyd.2", "DHX9-1", RNAseq_BRCA_isoforms$`Hybridization REF`)
RNAseq_BRCA_isoforms$`Hybridization REF` <- gsub("uc001gpr.2", "DHX9-2", RNAseq_BRCA_isoforms$`Hybridization REF`)


#make new data.frame with isoforms names as row.names
RNAseq_BRCA_isoforms <- data.frame(RNAseq_BRCA_isoforms, row.names = RNAseq_BRCA_isoforms$`Hybridization REF`)

#remove hybridization reference column
RNAseq_BRCA_isoforms$Hybridization.REF <- NULL

#transpose data.frame
RNAseq_BRCA_isoforms <- as.data.frame(t(RNAseq_BRCA_isoforms))

#make a new column for sample based on row.names
RNAseq_BRCA_isoforms$Sample <- row.names(RNAseq_BRCA_isoforms)

#split the sample identifier
out <- strsplit(as.character(RNAseq_BRCA_isoforms$Sample), "\\.")
test <- as.data.frame(do.call(rbind, out))

#build a patient_barcode
test$bcr_patient_barcode  <- do.call(paste, c(test[c("V1", "V2", "V3")], sep = "-")) 

#make a new column with the patient_barcode
RNAseq_BRCA_isoforms$patient_barcode <- test$bcr_patient_barcode

#build a sample barcode
test$Sample_Barcode <- do.call(paste, c(test[c("V1", "V2", "V3","V4")], sep = "-")) 

#make a new column wit the sample barcode
RNAseq_BRCA_isoforms$Sample_Barcode <- test$Sample_Barcode

#make a new column with sample type
RNAseq_BRCA_isoforms$Sample_Type <- test$V4

#rename the sample types as above
RNAseq_BRCA_isoforms$Sample_Type <- gsub("11A", "Normal", RNAseq_BRCA_isoforms$Sample_Type)
RNAseq_BRCA_isoforms$Sample_Type <- gsub("11B", "DROP", RNAseq_BRCA_isoforms$Sample_Type)

RNAseq_BRCA_isoforms$Sample_Type <- gsub("01A", "Primary", RNAseq_BRCA_isoforms$Sample_Type)
RNAseq_BRCA_isoforms$Sample_Type <- gsub("01B", "DROP", RNAseq_BRCA_isoforms$Sample_Type)

RNAseq_BRCA_isoforms$Sample_Type <- gsub("06A", "Metastatic", RNAseq_BRCA_isoforms$Sample_Type)
RNAseq_BRCA_isoforms$Sample_Type <- gsub("06B", "DROP", RNAseq_BRCA_isoforms$Sample_Type)

#remove B samples
RNAseq_BRCA_isoforms <- subset(RNAseq_BRCA_isoforms, !Sample_Type == "DROP")




#calculate total p150 expresssion, sum(p150-1,p150-2,p150-3)
RNAseq_BRCA_isoforms$p150 <- as.numeric(as.character(RNAseq_BRCA_isoforms$`p150-1`)) + 
  as.numeric(as.character(RNAseq_BRCA_isoforms$`p150-2`)) +
  as.numeric(as.character(RNAseq_BRCA_isoforms$`p150-3`))

RNAseq_BRCA_isoforms$DHX9 <- as.numeric(as.character(RNAseq_BRCA_isoforms$`DHX9-2`)) + 
  as.numeric(as.character(RNAseq_BRCA_isoforms$`DHX9-1`))


#make p110 numeric
RNAseq_BRCA_isoforms$p110 <- as.numeric(as.character(RNAseq_BRCA_isoforms$p110))

#subset only primary tumor samples
RNAseq_BRCA_isoforms_primary <- subset(RNAseq_BRCA_isoforms , Sample_Type == "Primary")

#merge with annotation data
RNAseq_BRCA_isoforms_primary <- merge(RNAseq_BRCA_isoforms_primary, Annotation, by.y = "BARCODE", by.x = "patient_barcode", all.x = TRUE)

#subset only normal samples
RNAseq_BRCA_isoforms_normal <- subset(RNAseq_BRCA_isoforms, Sample_Type == "Normal")

#fill in TNBC and PAM50 status for normal samples
RNAseq_BRCA_isoforms_normal$TNBC <- "Normal"
RNAseq_BRCA_isoforms_normal$PAM50 <- "Normal Breast"

#bind primary and normal
RNAseq_BRCA_isoforms <- rbind(RNAseq_BRCA_isoforms_normal, RNAseq_BRCA_isoforms_primary)

#plot correlation between p150 or p110 and DHX9 expression

ggplot(RNAseq_BRCA_isoforms, aes(DHX9, p150, colour = Sample_Type, alpha = Sample_Type)) + 
  geom_point() + theme_science() + scale_alpha_manual(values = c(0.5, 0.2), guide = "none") +
  geom_smooth(method = "lm") + labs(y = "ADAR p150 Expression\n (Normalized RSEM)", x = "DHX9 Expression\n (Normalized RSEM)", colour = "") +
  scale_color_manual(values = palette2, labels = c("Normal Breast", "Primary BRCA")) + 
  theme(legend.position = "bottom") + stat_cor(method="pearson", alpha = 1)
ggsave("DHX9_P150.tiff", height = 4.2, width = 4)
cor.test(RNAseq_BRCA_isoforms$DHX9, RNAseq_BRCA_isoforms$p150)

ggplot(RNAseq_BRCA_isoforms, aes(DHX9, p110, colour = Sample_Type, alpha = Sample_Type)) + 
  geom_point() + theme_science() + scale_alpha_manual(values = c(0.5, 0.2), guide = "none") +
  geom_smooth(method = "lm") + labs(y = "ADAR p110 Expression\n (Normalized RSEM)", x = "DHX9 Expression\n (Normalized RSEM)", colour = "") +
  scale_color_manual(values = palette2, labels = c("Normal Breast", "Primary BRCA")) + 
  theme(legend.position = "bottom") + stat_cor(method="pearson", alpha = 1)
ggsave("DHX9_P110.tiff", height = 4.2, width = 4)
cor.test(RNAseq_BRCA_isoforms$DHX9, RNAseq_BRCA_isoforms$p110)


