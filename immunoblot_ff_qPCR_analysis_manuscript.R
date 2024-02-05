library(ggplot2)
library(dplyr)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(data.table)
library(ggpubr)
library(viridis)
library(scales)
library(gmodels)
library(DescTools)
library(MetBrewer)
library(svglite)
library(ggbeeswarm)



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



scaleFUN <- function(x) sprintf("%.2f", x)

makeStars <- function(x){
  stars <- c("***", "**", "*", "ns")
  vec <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbviridis <- c("#440154FF", "#31688EFF", "#35B779FF", "#E69F00")

palette <- met.brewer(name="Veronese", n=7, type="discrete")

met.brewer(name="Veronese", n=7, type="discrete")

palette2 <- c(palette[7], palette[3])
pallete3 <- c(palette[1], palette[3], palette[5],palette[7])

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


setwd("/kacottre/Weber_Lab/APEX-ADAR_manuscript_dhx9/Source Data/")





#read western (immunoblot) data for knockdowns
immunoblot_skbr3_adar_ddx17 <- fread("immunoblot_skbr3_adar_ddx17.txt")

immunoblot_skbr3_adar_ddx17 <- na.omit(immunoblot_skbr3_adar_ddx17)

immunoblot_skbr3_adar_ddx17 <- subset(immunoblot_skbr3_adar_ddx17, !`fold_change` == "pending")


immunoblot_skbr3_adar_ddx17$sample <- paste(immunoblot_skbr3_adar_ddx17$shRNA_1, immunoblot_skbr3_adar_ddx17$shRNA_2, sep = "/")

immunoblot_skbr3_adar_ddx17$sample <- gsub("-", ".", immunoblot_skbr3_adar_ddx17$sample)

immunoblot_skbr3_adar_ddx17$fold_change <- as.numeric(immunoblot_skbr3_adar_ddx17$fold_change)

#group and summarise data
immunoblot_skbr3_adar_ddx17_sum <- dplyr::group_by(immunoblot_skbr3_adar_ddx17, sample, Gene)

immunoblot_skbr3_adar_ddx17_sum <- dplyr::summarise(immunoblot_skbr3_adar_ddx17_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_skbr3_adar_ddx17_sum$ci <- qnorm(0.975)*(immunoblot_skbr3_adar_ddx17_sum$sd_expression/sqrt(4))

#make plots and save

aov_ppkr <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx17, Gene == "p-PKR"))
tppkr <- TukeyHSD(aov_ppkr)
tppkr <- as.data.frame(tppkr$sample[,1:4])

out <- strsplit(as.character(row.names(tppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tppkr$stars <- makeStars(tppkr$`p adj`)

tppkr$xmin <- out$V1
tppkr$xmax <- out$V2
tppkr$y <- max(immunoblot_skbr3_adar_ddx17$fold_change[immunoblot_skbr3_adar_ddx17$Gene == "p-PKR"]*1.1)

tppkr <- tppkr[!grepl("shDDX17.1.*shDDX17.3|shDDX17.3.*shDDX17.1", rownames(tppkr)), ]

tppkr <- subset(tppkr, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx17, Gene == "p-PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX17.1", "shSCR/shDDX17.3", "shADAR/shSCR", "shADAR/shDDX17.1", "shADAR/shDDX17.3"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX17-1", "shSCR\nshDDX17-3", "shADAR1\nshSCR", "shADAR1\nshDDX17-1", "shADAR1\nshDDX17-3")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "p-PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "p-PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tppkr$xmin, 
               xmax = tppkr$xmax, 
               y.position = tppkr$y , label = tppkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_adar_ddx17_skbr3.tiff", height = 4.2, width = 4.2)




aov_peif2a <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx17, Gene == "p-eIF2a"))
tpeif2a <- TukeyHSD(aov_peif2a)
tpeif2a <- as.data.frame(tpeif2a$sample[,1:4])

out <- strsplit(as.character(row.names(tpeif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpeif2a$stars <- makeStars(tpeif2a$`p adj`)

tpeif2a$xmin <- out$V1
tpeif2a$xmax <- out$V2
tpeif2a$y <- max(immunoblot_skbr3_adar_ddx17$fold_change[immunoblot_skbr3_adar_ddx17$Gene == "p-eIF2a"]*1.1)

tpeif2a <- tpeif2a[!grepl("shDDX17.1.*shDDX17.3|shDDX17.3.*shDDX17.1", rownames(tpeif2a)), ]

tpeif2a <- subset(tpeif2a, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx17, Gene == "p-eIF2a"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX17.1", "shSCR/shDDX17.3", "shADAR/shSCR", "shADAR/shDDX17.1", "shADAR/shDDX17.3"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX17-1", "shSCR\nshDDX17-3", "shADAR1\nshSCR", "shADAR1\nshDDX17-1", "shADAR1\nshDDX17-3")) + 
  labs(x = "", y = "Relative p-eIF2a/eIF2a Abundance     ", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "p-eIF2a"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "p-eIF2a"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) 
ggsave("Western_plots/peif2a_adar_ddx17_skbr3.tiff", height = 4.2, width = 4.2)



aov_cparp <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx17, Gene == "c-PARP"))
tcparp <- TukeyHSD(aov_cparp)
tcparp <- as.data.frame(tcparp$sample[,1:4])

out <- strsplit(as.character(row.names(tcparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tcparp$stars <- makeStars(tcparp$`p adj`)

tcparp$xmin <- out$V1
tcparp$xmax <- out$V2
tcparp$y <- max(immunoblot_skbr3_adar_ddx17$fold_change[immunoblot_skbr3_adar_ddx17$Gene == "c-PARP"]*1.1)

tcparp <- tcparp[!grepl("shDDX17.1.*shDDX17.3|shDDX17.3.*shDDX17.1", rownames(tcparp)), ]

tcparp <- subset(tcparp, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx17, Gene == "c-PARP"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX17.1", "shSCR/shDDX17.3", "shADAR/shSCR", "shADAR/shDDX17.1", "shADAR/shDDX17.3"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX17-1", "shSCR\nshDDX17-3", "shADAR1\nshSCR", "shADAR1\nshDDX17-1", "shADAR1\nshDDX17-3")) + 
  labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "c-PARP"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "c-PARP"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tcparp$xmin, 
               xmax = tcparp$xmax, 
               y.position = tcparp$y , label = tcparp$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/cparp_adar_ddx17_skbr3.tiff", height = 4.2, width = 4.2)


aov_ddx17 <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx17, Gene == "DDX17"))
tddx17 <- TukeyHSD(aov_ddx17)
tddx17 <- as.data.frame(tddx17$sample[,1:4])

out <- strsplit(as.character(row.names(tddx17)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tddx17$stars <- makeStars(tddx17$`p adj`)

tddx17$xmin <- out$V1
tddx17$xmax <- out$V2
tddx17$y <- max(immunoblot_skbr3_adar_ddx17$fold_change[immunoblot_skbr3_adar_ddx17$Gene == "DDX17"]*1.1)

tddx17 <- tddx17[!grepl("shDDX17.1.*shDDX17.3|shDDX17.3.*shDDX17.1", rownames(tddx17)), ]

tddx17 <- subset(tddx17, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx17, Gene == "DDX17"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX17.1", "shSCR/shDDX17.3", "shADAR/shSCR", "shADAR/shDDX17.1", "shADAR/shDDX17.3"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX17-1", "shSCR\nshDDX17-3", "shADAR1\nshSCR", "shADAR1\nshDDX17-1", "shADAR1\nshDDX17-3")) + 
  labs(x = "", y = "Relative DDX17 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "DDX17"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "DDX17"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tddx17$xmin, 
               xmax = tddx17$xmax, 
               y.position = tddx17$y , label = tddx17$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ddx17_adar_ddx17_skbr3.tiff", height = 4.2, width = 4.2)


aov_adar <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx17, Gene == "ADARp110"))
tadar <- TukeyHSD(aov_adar)
tadar <- as.data.frame(tadar$sample[,1:4])

out <- strsplit(as.character(row.names(tadar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tadar$stars <- makeStars(tadar$`p adj`)

tadar$xmin <- out$V1
tadar$xmax <- out$V2
tadar$y <- max(immunoblot_skbr3_adar_ddx17$fold_change[immunoblot_skbr3_adar_ddx17$Gene == "ADARp110"]*1.1)

tadar <- tadar[!grepl("shDDX17.1.*shDDX17.3|shDDX17.3.*shDDX17.1", rownames(tadar)), ]

tadar <- subset(tadar, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx17, Gene == "ADARp110"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX17.1", "shSCR/shDDX17.3", "shADAR/shSCR", "shADAR/shDDX17.1", "shADAR/shDDX17.3"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX9-3", "shSCR\nshDDX9-5", "shADAR1\nshSCR", "shADAR1\nshDDX9-3", "shADAR1\nshDDX9-5")) + 
  labs(x = "", y = "Relative ADAR1-p110 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "ADARp110"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "ADARp110"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tadar$xmin, 
               xmax = tadar$xmax, 
               y.position = tadar$y , label = tadar$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/adar_adar_ddx17_skbr3.tiff", height = 4.2, width = 4.2)


aov_pkr <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_ddx17, Gene == "total PKR"))
tpkr <- TukeyHSD(aov_pkr)
tpkr <- as.data.frame(tpkr$sample[,1:4])

out <- strsplit(as.character(row.names(tpkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpkr$stars <- makeStars(tpkr$`p adj`)

tpkr$xmin <- out$V1
tpkr$xmax <- out$V2
tpkr$y <- max(immunoblot_skbr3_adar_ddx17$fold_change[immunoblot_skbr3_adar_ddx17$Gene == "total PKR"]*1.1)

tpkr <- tpkr[!grepl("shDDX17.1.*shDDX17.3|shDDX17.3.*shDDX17.1", rownames(tpkr)), ]

tpkr <- subset(tpkr, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_ddx17, Gene == "total PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX17.1", "shSCR/shDDX17.3", "shADAR/shSCR", "shADAR/shDDX17.1", "shADAR/shDDX17.3"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX9-3", "shSCR\nshDDX9-5", "shADAR1\nshSCR", "shADAR1\nshDDX9-3", "shADAR1\nshDDX9-5")) + 
  labs(x = "", y = "Relative total PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "total PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_ddx17_sum, Gene == "total PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) 
ggsave("Western_plots/pkr_adar_ddx17_skbr3.tiff", height = 4.2, width = 4.2)

aov_pkr <- NULL
tpkr <- NULL
aov_peif2a <- NULL
t_peif2a <- NULL
aov_ppkr <- NULL
tppkr <- NULL
aov_cparp <- NULL
tcparp <- NULL
aov_adar <- NULL
tadar <- NULL
aov_pkr <- NULL
tpkr <- NULL
aov_ddx17 <- NULL
tddx17 <- NULL

#read ff data for knockdowns
ff_skbr3_adar_ddx17 <- fread("DDX17 Foci Formation/ff_skbr3_ddx17_adar.txt")

ff_skbr3_adar_ddx17 <- na.omit(ff_skbr3_adar_ddx17)

ff_skbr3_adar_ddx17$sample <- paste(ff_skbr3_adar_ddx17$shRNA_1, ff_skbr3_adar_ddx17$shRNA_2, sep = "/")

ff_skbr3_adar_ddx17$sample <- gsub("-", ".", ff_skbr3_adar_ddx17$sample)

ff_skbr3_adar_ddx17$Relative_area <- as.numeric(ff_skbr3_adar_ddx17$Relative_area)

#group and summarise data
ff_skbr3_adar_ddx17_sum <- dplyr::group_by(ff_skbr3_adar_ddx17, sample)

ff_skbr3_adar_ddx17_sum <- dplyr::summarise(ff_skbr3_adar_ddx17_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_skbr3_adar_ddx17_sum$ci <- qnorm(0.975)*(ff_skbr3_adar_ddx17_sum$sd_area/sqrt(4))

aov_ff <- aov(Relative_area ~ sample, ff_skbr3_adar_ddx17)
tff <- TukeyHSD(aov_ff)
tff <- as.data.frame(tff$sample[,1:4])

out <- strsplit(as.character(row.names(tff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tff$stars <- makeStars(tff$`p adj`)

tff$xmin <- out$V1
tff$xmax <- out$V2
tff$y <- max(ff_skbr3_adar_ddx17$Relative_area*1.1)

tff <- tff[!grepl("shDDX17.1.*shDDX17.3|shDDX17.3.*shDDX17.1", rownames(tff)), ]

tff <- subset(tff, !stars == "ns")

ggplot(ff_skbr3_adar_ddx17, aes(sample, Relative_area)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDDX17.1", "shSCR/shDDX17.3", "shADAR/shSCR", "shADAR/shDDX17.1", "shADAR/shDDX17.3"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDDX9-3", "shSCR\nshDDX9-5", "shADAR1\nshSCR", "shADAR1\nshDDX9-3", "shADAR1\nshDDX9-5")) + 
  labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = ff_skbr3_adar_ddx17_sum, aes(sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = ff_skbr3_adar_ddx17_sum, aes(x = sample, y = mean_area, ymin = mean_area - sd_area, ymax = mean_area + sd_area), width=0.2) 
ggsave("FF_plots/adar_ddx17_skbr3.tiff", height = 4.2, width = 4.2, units = "in")


aov_ff <- NULL
tff <- NULL


#read western (immunoblot) data for knockdowns
immunoblot_skbr3_adar_dhx9 <- fread("immunoblot_skbr3_adar_dhx9.txt")

immunoblot_skbr3_adar_dhx9 <- na.omit(immunoblot_skbr3_adar_dhx9)

immunoblot_skbr3_adar_dhx9 <- subset(immunoblot_skbr3_adar_dhx9, !`fold_change` == "Pending")


immunoblot_skbr3_adar_dhx9$sample <- paste(immunoblot_skbr3_adar_dhx9$shRNA_1, immunoblot_skbr3_adar_dhx9$shRNA_2, sep = "/")

immunoblot_skbr3_adar_dhx9$sample <- gsub("-", ".", immunoblot_skbr3_adar_dhx9$sample)

immunoblot_skbr3_adar_dhx9$fold_change <- as.numeric(immunoblot_skbr3_adar_dhx9$fold_change)

#group and summarise data
immunoblot_skbr3_adar_dhx9_sum <- dplyr::group_by(immunoblot_skbr3_adar_dhx9, sample, Gene)

immunoblot_skbr3_adar_dhx9_sum <- dplyr::summarise(immunoblot_skbr3_adar_dhx9_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_skbr3_adar_dhx9_sum$ci <- qnorm(0.975)*(immunoblot_skbr3_adar_dhx9_sum$sd_expression/sqrt(4))

#make plots and save

aov_ppkr <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_dhx9, Gene == "p-PKR"))
tppkr <- TukeyHSD(aov_ppkr)
tppkr <- as.data.frame(tppkr$sample[,1:4])

out <- strsplit(as.character(row.names(tppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tppkr$stars <- makeStars(tppkr$`p adj`)

tppkr$xmin <- out$V1
tppkr$xmax <- out$V2
tppkr$y <- max(immunoblot_skbr3_adar_dhx9$fold_change[immunoblot_skbr3_adar_dhx9$Gene == "p-PKR"]*1.1)

tppkr <- tppkr[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tppkr)), ]

tppkr <- subset(tppkr, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_dhx9, Gene == "p-PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "p-PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "p-PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tppkr$xmin, 
               xmax = tppkr$xmax, 
               y.position = tppkr$y , label = tppkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_adar_dhx9_skbr3.tiff", height = 4.2, width = 4.2)


aov_peif2a <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_dhx9, Gene == "p-eIF2a"))
tpeif2a <- TukeyHSD(aov_peif2a)
tpeif2a <- as.data.frame(tpeif2a$sample[,1:4])

out <- strsplit(as.character(row.names(tpeif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpeif2a$stars <- makeStars(tpeif2a$`p adj`)

tpeif2a$xmin <- out$V1
tpeif2a$xmax <- out$V2
tpeif2a$y <- max(immunoblot_skbr3_adar_dhx9$fold_change[immunoblot_skbr3_adar_dhx9$Gene == "p-eIF2a"]*1.1)

tpeif2a <- tpeif2a[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tpeif2a)), ]

tpeif2a <- subset(tpeif2a, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_dhx9, Gene == "p-eIF2a"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative p-eIF2a/eIF2a Abundance     ", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "p-eIF2a"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "p-eIF2a"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tpeif2a$xmin, 
               xmax = tpeif2a$xmax, 
               y.position = tpeif2a$y , label = tpeif2a$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/peif2a_adar_dhx9_skbr3.tiff", height = 4.2, width = 4.2)


aov_cparp <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_dhx9, Gene == "c-PARP"))
tcparp <- TukeyHSD(aov_cparp)
tcparp <- as.data.frame(tcparp$sample[,1:4])

out <- strsplit(as.character(row.names(tcparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tcparp$stars <- makeStars(tcparp$`p adj`)

tcparp$xmin <- out$V1
tcparp$xmax <- out$V2
tcparp$y <- max(immunoblot_skbr3_adar_dhx9$fold_change[immunoblot_skbr3_adar_dhx9$Gene == "c-PARP"]*1.1)

tcparp <- tcparp[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tcparp)), ]

tcparp <- subset(tcparp, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_dhx9, Gene == "c-PARP"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "c-PARP"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "c-PARP"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tcparp$xmin, 
               xmax = tcparp$xmax, 
               y.position = tcparp$y , label = tcparp$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/cparp_adar_dhx9_skbr3.tiff", height = 4.2, width = 4.2)


aov_dhx9 <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_dhx9, Gene == "DHX9"))
tdhx9 <- TukeyHSD(aov_dhx9)
tdhx9 <- as.data.frame(tdhx9$sample[,1:4])

out <- strsplit(as.character(row.names(tdhx9)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tdhx9$stars <- makeStars(tdhx9$`p adj`)

tdhx9$xmin <- out$V1
tdhx9$xmax <- out$V2
tdhx9$y <- max(immunoblot_skbr3_adar_dhx9$fold_change[immunoblot_skbr3_adar_dhx9$Gene == "DHX9"]*1.1)

tdhx9 <- tdhx9[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tdhx9)), ]

tdhx9 <- subset(tdhx9, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_dhx9, Gene == "DHX9"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative DHX9 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "DHX9"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "DHX9"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tdhx9$xmin, 
               xmax = tdhx9$xmax, 
               y.position = tdhx9$y , label = tdhx9$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/dhx9_adar_dhx9_skbr3.tiff", height = 4.2, width = 4.2)




aov_adar <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_dhx9, Gene == "ADARp110"))
tadar <- TukeyHSD(aov_adar)
tadar <- as.data.frame(tadar$sample[,1:4])

out <- strsplit(as.character(row.names(tadar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tadar$stars <- makeStars(tadar$`p adj`)

tadar$xmin <- out$V1
tadar$xmax <- out$V2
tadar$y <- max(immunoblot_skbr3_adar_dhx9$fold_change[immunoblot_skbr3_adar_dhx9$Gene == "ADARp110"]*1.1)

tadar <- tadar[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tadar)), ]

tadar <- subset(tadar, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_dhx9, Gene == "ADARp110"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative ADAR1-p110 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "ADARp110"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "ADARp110"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tadar$xmin, 
               xmax = tadar$xmax, 
               y.position = tadar$y , label = tadar$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/adar_adar_dhx9_skbr3.tiff", height = 4.2, width = 4.2)


aov_pkr <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_dhx9, Gene == "total PKR"))
tpkr <- TukeyHSD(aov_pkr)
tpkr <- as.data.frame(tpkr$sample[,1:4])

out <- strsplit(as.character(row.names(tpkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpkr$stars <- makeStars(tpkr$`p adj`)

tpkr$xmin <- out$V1
tpkr$xmax <- out$V2
tpkr$y <- max(immunoblot_skbr3_adar_dhx9$fold_change[immunoblot_skbr3_adar_dhx9$Gene == "total PKR"]*1.1)

tpkr <- tpkr[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tpkr)), ]

tpkr <- subset(tpkr, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_dhx9, Gene == "total PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative total PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "total PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "total PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tpkr$xmin, 
               xmax = tpkr$xmax, 
               y.position = tpkr$y , label = tpkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/pkr_adar_dhx9_skbr3.tiff", height = 4.2, width = 4.2)




aov_ISG15 <- aov(fold_change ~ sample, subset(immunoblot_skbr3_adar_dhx9, Gene == "ISG15"))
tISG15 <- TukeyHSD(aov_ISG15)
tISG15 <- as.data.frame(tISG15$sample[,1:4])

out <- strsplit(as.character(row.names(tISG15)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tISG15$stars <- makeStars(tISG15$`p adj`)

tISG15$xmin <- out$V1
tISG15$xmax <- out$V2
tISG15$y <- max(immunoblot_skbr3_adar_dhx9$fold_change[immunoblot_skbr3_adar_dhx9$Gene == "ISG15"]*1.1)

tISG15 <- tISG15[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tISG15)), ]

tISG15 <- subset(tISG15, !stars == "ns")

ggplot(subset(immunoblot_skbr3_adar_dhx9, Gene == "ISG15"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative ISG15 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "ISG15"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_adar_dhx9_sum, Gene == "ISG15"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tISG15$xmin, 
               xmax = tISG15$xmax, 
               y.position = tISG15$y , label = tISG15$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ISG15_adar_dhx9_skbr3.tiff", height = 4.2, width = 4.2)


aov_pkr <- NULL
tpkr <- NULL
aov_peif2a <- NULL
t_peif2a <- NULL
aov_ppkr <- NULL
tppkr <- NULL
aov_cparp <- NULL
tcparp <- NULL
aov_adar <- NULL
tadar <- NULL
aov_pkr <- NULL
tpkr <- NULL
aov_dhx9 <- NULL
tdhx9 <- NULL
aov_ISG15 <- NULL
tISG15 <- NULL


#read ff data for knockdowns
ff_skbr3_adar_dhx9 <- fread("DHX9 Foci Formation/SKRB3 ADAR DHX9/ff_skbr3_dhx9_adar.txt")

ff_skbr3_adar_dhx9 <- na.omit(ff_skbr3_adar_dhx9)

ff_skbr3_adar_dhx9$sample <- paste(ff_skbr3_adar_dhx9$shRNA_1, ff_skbr3_adar_dhx9$shRNA_2, sep = "/")

ff_skbr3_adar_dhx9$sample <- gsub("-", ".", ff_skbr3_adar_dhx9$sample)

ff_skbr3_adar_dhx9$Relative_area <- as.numeric(ff_skbr3_adar_dhx9$Relative_area)

#group and summarise data
ff_skbr3_adar_dhx9_sum <- dplyr::group_by(ff_skbr3_adar_dhx9, sample)

ff_skbr3_adar_dhx9_sum <- dplyr::summarise(ff_skbr3_adar_dhx9_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_skbr3_adar_dhx9_sum$ci <- qnorm(0.975)*(ff_skbr3_adar_dhx9_sum$sd_area/sqrt(4))




aov_ff <- aov(Relative_area ~ sample, ff_skbr3_adar_dhx9)
tff <- TukeyHSD(aov_ff)
tff <- as.data.frame(tff$sample[,1:4])

out <- strsplit(as.character(row.names(tff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tff$stars <- makeStars(tff$`p adj`)

tff$xmin <- out$V1
tff$xmax <- out$V2
tff$y <- max(ff_skbr3_adar_dhx9$Relative_area*1.1)

tff <- tff[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tff)), ]

tff <- subset(tff, !stars == "ns")

ggplot(ff_skbr3_adar_dhx9, aes(sample, Relative_area)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = ff_skbr3_adar_dhx9_sum, aes(sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = ff_skbr3_adar_dhx9_sum, aes(x = sample, y = mean_area, ymin = mean_area - sd_area, ymax = mean_area + sd_area), width=0.2) +
  geom_bracket(xmin = tff$xmin, 
               xmax = tff$xmax, 
               y.position = tff$y , label = tff$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("FF_plots/adar_dhx9_skbr3.tiff", height = 4.2, width = 4.2, units = "in")


aov_ff <- NULL
tff <- NULL



#read western (immunoblot) data for knockdowns
immunoblot_mcf7_adar_dhx9 <- fread("immunoblot_mcf7_adar_dhx9.txt")

immunoblot_mcf7_adar_dhx9 <- na.omit(immunoblot_mcf7_adar_dhx9)

immunoblot_mcf7_adar_dhx9 <- subset(immunoblot_mcf7_adar_dhx9, !`fold_change` == "Pending")


immunoblot_mcf7_adar_dhx9$sample <- paste(immunoblot_mcf7_adar_dhx9$shRNA_1, immunoblot_mcf7_adar_dhx9$shRNA_2, sep = "/")

immunoblot_mcf7_adar_dhx9$sample <- gsub("-", ".", immunoblot_mcf7_adar_dhx9$sample)

immunoblot_mcf7_adar_dhx9$fold_change <- as.numeric(immunoblot_mcf7_adar_dhx9$fold_change)

#group and summarise data
immunoblot_mcf7_adar_dhx9_sum <- dplyr::group_by(immunoblot_mcf7_adar_dhx9, sample, Gene)

immunoblot_mcf7_adar_dhx9_sum <- dplyr::summarise(immunoblot_mcf7_adar_dhx9_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_mcf7_adar_dhx9_sum$ci <- qnorm(0.975)*(immunoblot_mcf7_adar_dhx9_sum$sd_expression/sqrt(4))

#make plots and save

aov_ppkr <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_dhx9, Gene == "p-PKR"))
tppkr <- TukeyHSD(aov_ppkr)
tppkr <- as.data.frame(tppkr$sample[,1:4])

out <- strsplit(as.character(row.names(tppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tppkr$stars <- makeStars(tppkr$`p adj`)

tppkr$xmin <- out$V1
tppkr$xmax <- out$V2
tppkr$y <- max(immunoblot_mcf7_adar_dhx9$fold_change[immunoblot_mcf7_adar_dhx9$Gene == "p-PKR"]*1.1)

tppkr <- tppkr[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tppkr)), ]

tppkr <- subset(tppkr, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_dhx9, Gene == "p-PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "p-PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "p-PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tppkr$xmin, 
               xmax = tppkr$xmax, 
               y.position = tppkr$y , label = tppkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_adar_dhx9_mcf7.tiff", height = 4.2, width = 4.2)


aov_peif2a <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_dhx9, Gene == "p-eIF2a"))
tpeif2a <- TukeyHSD(aov_peif2a)
tpeif2a <- as.data.frame(tpeif2a$sample[,1:4])

out <- strsplit(as.character(row.names(tpeif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpeif2a$stars <- makeStars(tpeif2a$`p adj`)

tpeif2a$xmin <- out$V1
tpeif2a$xmax <- out$V2
tpeif2a$y <- max(immunoblot_mcf7_adar_dhx9$fold_change[immunoblot_mcf7_adar_dhx9$Gene == "p-eIF2a"]*1.1)

tpeif2a <- tpeif2a[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tpeif2a)), ]

tpeif2a <- subset(tpeif2a, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_dhx9, Gene == "p-eIF2a"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative p-eIF2a/eIF2a Abundance     ", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "p-eIF2a"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "p-eIF2a"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tpeif2a$xmin, 
               xmax = tpeif2a$xmax, 
               y.position = tpeif2a$y , label = tpeif2a$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/peif2a_adar_dhx9_mcf7.tiff", height = 4.2, width = 4.2)


aov_cparp <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_dhx9, Gene == "c-PARP"))
tcparp <- TukeyHSD(aov_cparp)
tcparp <- as.data.frame(tcparp$sample[,1:4])

out <- strsplit(as.character(row.names(tcparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tcparp$stars <- makeStars(tcparp$`p adj`)

tcparp$xmin <- out$V1
tcparp$xmax <- out$V2
tcparp$y <- max(immunoblot_mcf7_adar_dhx9$fold_change[immunoblot_mcf7_adar_dhx9$Gene == "c-PARP"]*1.1)

tcparp <- tcparp[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tcparp)), ]

tcparp <- subset(tcparp, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_dhx9, Gene == "c-PARP"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "c-PARP"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "c-PARP"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2)# +
  #geom_bracket(xmin = tcparp$xmin, 
            #   xmax = tcparp$xmax, 
           #    y.position = tcparp$y , label = tcparp$stars, 
             #  tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/cparp_adar_dhx9_mcf7.tiff", height = 4.2, width = 4.2)


aov_dhx9 <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_dhx9, Gene == "DHX9"))
tdhx9 <- TukeyHSD(aov_dhx9)
tdhx9 <- as.data.frame(tdhx9$sample[,1:4])

out <- strsplit(as.character(row.names(tdhx9)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tdhx9$stars <- makeStars(tdhx9$`p adj`)

tdhx9$xmin <- out$V1
tdhx9$xmax <- out$V2
tdhx9$y <- max(immunoblot_mcf7_adar_dhx9$fold_change[immunoblot_mcf7_adar_dhx9$Gene == "DHX9"]*1.1)

tdhx9 <- tdhx9[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tdhx9)), ]

tdhx9 <- subset(tdhx9, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_dhx9, Gene == "DHX9"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative DHX9 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "DHX9"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "DHX9"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tdhx9$xmin, 
               xmax = tdhx9$xmax, 
               y.position = tdhx9$y , label = tdhx9$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/dhx9_adar_dhx9_mcf7.tiff", height = 4.2, width = 4.2)




aov_adar <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_dhx9, Gene == "ADARp110"))
tadar <- TukeyHSD(aov_adar)
tadar <- as.data.frame(tadar$sample[,1:4])

out <- strsplit(as.character(row.names(tadar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tadar$stars <- makeStars(tadar$`p adj`)

tadar$xmin <- out$V1
tadar$xmax <- out$V2
tadar$y <- max(immunoblot_mcf7_adar_dhx9$fold_change[immunoblot_mcf7_adar_dhx9$Gene == "ADARp110"]*1.1)

tadar <- tadar[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tadar)), ]

tadar <- subset(tadar, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_dhx9, Gene == "ADARp110"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative ADAR1-p110 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "ADARp110"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "ADARp110"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tadar$xmin, 
               xmax = tadar$xmax, 
               y.position = tadar$y , label = tadar$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/adar_adar_dhx9_mcf7.tiff", height = 4.2, width = 4.2)


aov_pkr <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_dhx9, Gene == "total PKR"))
tpkr <- TukeyHSD(aov_pkr)
tpkr <- as.data.frame(tpkr$sample[,1:4])

out <- strsplit(as.character(row.names(tpkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tpkr$stars <- makeStars(tpkr$`p adj`)

tpkr$xmin <- out$V1
tpkr$xmax <- out$V2
tpkr$y <- max(immunoblot_mcf7_adar_dhx9$fold_change[immunoblot_mcf7_adar_dhx9$Gene == "total PKR"]*1.1)

tpkr <- tpkr[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tpkr)), ]

tpkr <- subset(tpkr, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_dhx9, Gene == "total PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative total PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "total PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "total PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tpkr$xmin, 
               xmax = tpkr$xmax, 
               y.position = tpkr$y , label = tpkr$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/pkr_adar_dhx9_mcf7.tiff", height = 4.2, width = 4.2)



aov_ISG15 <- aov(fold_change ~ sample, subset(immunoblot_mcf7_adar_dhx9, Gene == "ISG15"))
tISG15 <- TukeyHSD(aov_ISG15)
tISG15 <- as.data.frame(tISG15$sample[,1:4])

out <- strsplit(as.character(row.names(tISG15)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tISG15$stars <- makeStars(tISG15$`p adj`)

tISG15$xmin <- out$V1
tISG15$xmax <- out$V2
tISG15$y <- max(immunoblot_mcf7_adar_dhx9$fold_change[immunoblot_mcf7_adar_dhx9$Gene == "ISG15"]*1.1)

tISG15 <- tISG15[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tISG15)), ]

tISG15 <- subset(tISG15, !stars == "ns")

ggplot(subset(immunoblot_mcf7_adar_dhx9, Gene == "ISG15"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative ISG15 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "ISG15"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_adar_dhx9_sum, Gene == "ISG15"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = tISG15$xmin, 
               xmax = tISG15$xmax, 
               y.position = tISG15$y , label = tISG15$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ISG15_adar_dhx9_mcf7.tiff", height = 4.2, width = 4.2)

aov_pkr <- NULL
tpkr <- NULL
aov_peif2a <- NULL
t_peif2a <- NULL
aov_ppkr <- NULL
tppkr <- NULL
aov_cparp <- NULL
tcparp <- NULL
aov_adar <- NULL
tadar <- NULL
aov_pkr <- NULL
tpkr <- NULL
aov_dhx9 <- NULL
tdhx9 <- NULL
aov_ISG15 <- NULL
tISG15 <- NULL


#read ff data for knockdowns
ff_mcf7_adar_dhx9 <- fread("DHX9 Foci Formation/MCF7/ff_mcf7_dhx9_adar.txt")

ff_mcf7_adar_dhx9 <- na.omit(ff_mcf7_adar_dhx9)

ff_mcf7_adar_dhx9$sample <- paste(ff_mcf7_adar_dhx9$shRNA_1, ff_mcf7_adar_dhx9$shRNA_2, sep = "/")

ff_mcf7_adar_dhx9$sample <- gsub("-", ".", ff_mcf7_adar_dhx9$sample)

ff_mcf7_adar_dhx9$Relative_area <- as.numeric(ff_mcf7_adar_dhx9$Relative_area)*100

#group and summarise data
ff_mcf7_adar_dhx9_sum <- dplyr::group_by(ff_mcf7_adar_dhx9, sample)

ff_mcf7_adar_dhx9_sum <- dplyr::summarise(ff_mcf7_adar_dhx9_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_mcf7_adar_dhx9_sum$ci <- qnorm(0.975)*(ff_mcf7_adar_dhx9_sum$sd_area/sqrt(4))




aov_ff <- aov(Relative_area ~ sample, ff_mcf7_adar_dhx9)
tff <- TukeyHSD(aov_ff)
tff <- as.data.frame(tff$sample[,1:4])

out <- strsplit(as.character(row.names(tff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

tff$stars <- makeStars(tff$`p adj`)

tff$xmin <- out$V1
tff$xmax <- out$V2
tff$y <- max(ff_mcf7_adar_dhx9$Relative_area*1.1)

tff <- tff[!grepl("shDHX9.3.*shDHX9.5|shDHX9.5.*shDHX9.3", rownames(tff)), ]

tff <- subset(tff, !stars == "ns")

ggplot(ff_mcf7_adar_dhx9, aes(sample, Relative_area)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/shSCR", "shSCR/shDHX9.3", "shSCR/shDHX9.5", "shADAR/shSCR", "shADAR/shDHX9.3", "shADAR/shDHX9.5"),
                   labels = c("shSCR\nshSCR", "shSCR\nshDHX9-3", "shSCR\nshDHX9-5", "shADAR1\nshSCR", "shADAR1\nshDHX9-3", "shADAR1\nshDHX9-5")) + 
  labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = ff_mcf7_adar_dhx9_sum, aes(sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = ff_mcf7_adar_dhx9_sum, aes(x = sample, y = mean_area, ymin = mean_area - sd_area, ymax = mean_area + sd_area), width=0.2) +
  geom_bracket(xmin = tff$xmin, 
               xmax = tff$xmax, 
               y.position = tff$y , label = tff$stars, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("FF_plots/adar_dhx9_mcf7.tiff", height = 4.2, width = 4.2, units = "in")


aov_ff <- NULL
tff <- NULL




####Immunoblot TNBC lines
immunoblot_tnbc_dhx9 <- fread("immunoblot_TNBC_dhx9.txt")

immunoblot_tnbc_dhx9 <- na.omit(immunoblot_tnbc_dhx9)

immunoblot_tnbc_dhx9$shRNA <- gsub("-", ".", immunoblot_tnbc_dhx9$shRNA)

#group and summarise data
immunoblot_tnbc_dhx9_sum <- dplyr::group_by(immunoblot_tnbc_dhx9, Cell_line, Gene, shRNA)

immunoblot_tnbc_dhx9_sum<- dplyr::summarise(immunoblot_tnbc_dhx9_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

tnbc_dhx9_ppkr_mb231 <-subset(immunoblot_tnbc_dhx9, Gene == "p-PKR" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_ppkr_mb231)
summary(mb231)

tnbc_dhx9_ppkr_mb231$shRNA <- factor(as.factor(tnbc_dhx9_ppkr_mb231$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb231ppkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_ppkr_mb231)
mb231ppkr <- as.data.frame(mb231ppkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231ppkr$group1 <- out$V1
mb231ppkr$group2 <- out$V2
mb231ppkr$pval <- makeStars(mb231ppkr$pval)
mb231ppkr$y.position <- max(tnbc_dhx9_ppkr_mb231$fold_change*1.1)
mb231ppkr$Cell_line = "MB231"

tnbc_dhx9_ppkr_hcc1806 <-subset(immunoblot_tnbc_dhx9, Gene == "p-PKR" & Cell_line == "HCC1806")
hcc1806 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_ppkr_hcc1806)
summary(hcc1806)

tnbc_dhx9_ppkr_hcc1806$shRNA <- factor(as.factor(tnbc_dhx9_ppkr_hcc1806$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

hcc1806ppkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_ppkr_hcc1806)
hcc1806ppkr <- as.data.frame(hcc1806ppkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(hcc1806ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806ppkr$group1 <- out$V1
hcc1806ppkr$group2 <- out$V2
hcc1806ppkr$pval <- makeStars(hcc1806ppkr$pval)
hcc1806ppkr$y.position <- max(tnbc_dhx9_ppkr_hcc1806$fold_change*1.1)
hcc1806ppkr$Cell_line = "HCC1806"

tnbc_dhx9_ppkr_bt549 <-subset(immunoblot_tnbc_dhx9, Gene == "p-PKR" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_ppkr_bt549)
summary(bt549)

tnbc_dhx9_ppkr_bt549$shRNA <- factor(as.factor(tnbc_dhx9_ppkr_bt549$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

bt549ppkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_ppkr_bt549)
bt549ppkr <- as.data.frame(bt549ppkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ppkr$group1 <- out$V1
bt549ppkr$group2 <- out$V2
bt549ppkr$pval <- makeStars(bt549ppkr$pval)
bt549ppkr$y.position <- max(tnbc_dhx9_ppkr_bt549$fold_change*1.1)
bt549ppkr$Cell_line = "BT549"


tnbc_dhx9_ppkr_mb468 <-subset(immunoblot_tnbc_dhx9, Gene == "p-PKR" & Cell_line == "MB468")
mb468 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_ppkr_mb468)
summary(mb468)

tnbc_dhx9_ppkr_mb468$shRNA <- factor(as.factor(tnbc_dhx9_ppkr_mb468$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb468ppkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_ppkr_mb468)
mb468ppkr <- as.data.frame(mb468ppkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb468ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468ppkr$group1 <- out$V1
mb468ppkr$group2 <- out$V2
mb468ppkr$pval <- makeStars(mb468ppkr$pval)
mb468ppkr$y.position <- max(tnbc_dhx9_ppkr_mb468$fold_change*1.1)
mb468ppkr$Cell_line = "MB468"

tnbcppkr <- rbind(mb468ppkr, mb231ppkr, hcc1806ppkr, bt549ppkr)
tnbcppkr <- subset(tnbcppkr, !pval == "ns")

ggplot(subset(immunoblot_tnbc_dhx9, Gene == "p-PKR"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDHX9.3", "shDHX9.5"), labels = c("shSCR", "shDHX9-3", "shDHX9-5")) + labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + 
  scale_colour_manual(values = pallete3, 
                      labels = c(BT549 = "BT-549", HCC1806  = "HCC1806", MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")) +
  geom_col(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "p-PKR", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "p-PKR", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcppkr, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/ppkr_tnbc_dhx9_facet.tiff", height =5, width = 5, units = "in")


tnbc_dhx9_cparp_mb231 <-subset(immunoblot_tnbc_dhx9, Gene == "c-PARP" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_cparp_mb231)
summary(mb231)

tnbc_dhx9_cparp_mb231$shRNA <- factor(as.factor(tnbc_dhx9_cparp_mb231$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb231cparp <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_cparp_mb231)
mb231cparp <- as.data.frame(mb231cparp$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231cparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231cparp$group1 <- out$V1
mb231cparp$group2 <- out$V2
mb231cparp$pval <- makeStars(mb231cparp$pval)
mb231cparp$y.position <- max(tnbc_dhx9_cparp_mb231$fold_change*1.1)
mb231cparp$Cell_line = "MB231"

tnbc_dhx9_cparp_hcc1806 <-subset(immunoblot_tnbc_dhx9, Gene == "c-PARP" & Cell_line == "HCC1806")
hcc1806 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_cparp_hcc1806)
summary(hcc1806)

tnbc_dhx9_cparp_hcc1806$shRNA <- factor(as.factor(tnbc_dhx9_cparp_hcc1806$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

hcc1806cparp <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_cparp_hcc1806)
hcc1806cparp <- as.data.frame(hcc1806cparp$shSCR[,1:4])

out <- strsplit(as.character(row.names(hcc1806cparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806cparp$group1 <- out$V1
hcc1806cparp$group2 <- out$V2
hcc1806cparp$pval <- makeStars(hcc1806cparp$pval)
hcc1806cparp$y.position <- max(tnbc_dhx9_cparp_hcc1806$fold_change*1.1)
hcc1806cparp$Cell_line = "HCC1806"

tnbc_dhx9_cparp_bt549 <-subset(immunoblot_tnbc_dhx9, Gene == "c-PARP" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_cparp_bt549)
summary(bt549)

tnbc_dhx9_cparp_bt549$shRNA <- factor(as.factor(tnbc_dhx9_cparp_bt549$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

bt549cparp <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_cparp_bt549)
bt549cparp <- as.data.frame(bt549cparp$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549cparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549cparp$group1 <- out$V1
bt549cparp$group2 <- out$V2
bt549cparp$pval <- makeStars(bt549cparp$pval)
bt549cparp$y.position <- max(tnbc_dhx9_cparp_bt549$fold_change*1.1)
bt549cparp$Cell_line = "BT549"


tnbc_dhx9_cparp_mb468 <-subset(immunoblot_tnbc_dhx9, Gene == "c-PARP" & Cell_line == "MB468")
mb468 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_cparp_mb468)
summary(mb468)

tnbc_dhx9_cparp_mb468$shRNA <- factor(as.factor(tnbc_dhx9_cparp_mb468$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb468cparp <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_cparp_mb468)
mb468cparp <- as.data.frame(mb468cparp$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb468cparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468cparp$group1 <- out$V1
mb468cparp$group2 <- out$V2
mb468cparp$pval <- makeStars(mb468cparp$pval)
mb468cparp$y.position <- max(tnbc_dhx9_cparp_mb468$fold_change*1.1)
mb468cparp$Cell_line = "MB468"

tnbccparp <- rbind(mb468cparp, mb231cparp, hcc1806cparp, bt549cparp)
tnbccparp <- subset(tnbccparp, !pval == "ns")

ggplot(subset(immunoblot_tnbc_dhx9, Gene == "c-PARP"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDHX9.3", "shDHX9.5"), labels = c("shSCR", "shDHX9-3", "shDHX9-5")) + labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + scale_colour_manual(values = pallete3, 
                                                                                    labels = c(BT549 = "BT-549", HCC1806  = "HCC1806", MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")) +
  geom_col(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "c-PARP", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "c-PARP", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbccparp, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/cparp_tnbc_dhx9_facet.tiff", height =5, width = 5, units = "in")


tnbc_dhx9_peif2a_mb231 <-subset(immunoblot_tnbc_dhx9, Gene == "p-eIF2a" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_peif2a_mb231)
summary(mb231)

tnbc_dhx9_peif2a_mb231$shRNA <- factor(as.factor(tnbc_dhx9_peif2a_mb231$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb231peif2a <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_peif2a_mb231)
mb231peif2a <- as.data.frame(mb231peif2a$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231peif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231peif2a$group1 <- out$V1
mb231peif2a$group2 <- out$V2
mb231peif2a$pval <- makeStars(mb231peif2a$pval)
mb231peif2a$y.position <- max(tnbc_dhx9_peif2a_mb231$fold_change*1.1)
mb231peif2a$Cell_line = "MB231"

tnbc_dhx9_peif2a_hcc1806 <-subset(immunoblot_tnbc_dhx9, Gene == "p-eIF2a" & Cell_line == "HCC1806")
hcc1806 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_peif2a_hcc1806)
summary(hcc1806)

tnbc_dhx9_peif2a_hcc1806$shRNA <- factor(as.factor(tnbc_dhx9_peif2a_hcc1806$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

hcc1806peif2a <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_peif2a_hcc1806)
hcc1806peif2a <- as.data.frame(hcc1806peif2a$shSCR[,1:4])

out <- strsplit(as.character(row.names(hcc1806peif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806peif2a$group1 <- out$V1
hcc1806peif2a$group2 <- out$V2
hcc1806peif2a$pval <- makeStars(hcc1806peif2a$pval)
hcc1806peif2a$y.position <- max(tnbc_dhx9_peif2a_hcc1806$fold_change*1.1)
hcc1806peif2a$Cell_line = "HCC1806"

tnbc_dhx9_peif2a_bt549 <-subset(immunoblot_tnbc_dhx9, Gene == "p-eIF2a" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_peif2a_bt549)
summary(bt549)

tnbc_dhx9_peif2a_bt549$shRNA <- factor(as.factor(tnbc_dhx9_peif2a_bt549$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

bt549peif2a <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_peif2a_bt549)
bt549peif2a <- as.data.frame(bt549peif2a$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549peif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549peif2a$group1 <- out$V1
bt549peif2a$group2 <- out$V2
bt549peif2a$pval <- makeStars(bt549peif2a$pval)
bt549peif2a$y.position <- max(tnbc_dhx9_peif2a_bt549$fold_change*1.1)
bt549peif2a$Cell_line = "BT549"


tnbc_dhx9_peif2a_mb468 <-subset(immunoblot_tnbc_dhx9, Gene == "p-eIF2a" & Cell_line == "MB468")
mb468 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_peif2a_mb468)
summary(mb468)

tnbc_dhx9_peif2a_mb468$shRNA <- factor(as.factor(tnbc_dhx9_peif2a_mb468$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb468peif2a <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_peif2a_mb468)
mb468peif2a <- as.data.frame(mb468peif2a$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb468peif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468peif2a$group1 <- out$V1
mb468peif2a$group2 <- out$V2
mb468peif2a$pval <- makeStars(mb468peif2a$pval)
mb468peif2a$y.position <- max(tnbc_dhx9_peif2a_mb468$fold_change*1.1)
mb468peif2a$Cell_line = "MB468"

tnbcpeif2a <- rbind(mb468peif2a, mb231peif2a, hcc1806peif2a, bt549peif2a)
tnbcpeif2a <- subset(tnbcpeif2a, !pval == "ns")

ggplot(subset(immunoblot_tnbc_dhx9, Gene == "p-eIF2a"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDHX9.3", "shDHX9.5"), labels = c("shSCR", "shDHX9-3", "shDHX9-5")) + 
  labs(x = "", y = "Relative p-eIF2a/eIF2a Abundance     ", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + scale_colour_manual(values = pallete3, 
                                                                                    labels = c(BT549 = "BT-549", HCC1806  = "HCC1806", MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")) +
  geom_col(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "p-eIF2a", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "p-eIF2a", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcpeif2a, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/peif2a_tnbc_dhx9_facet.tiff", height =5, width = 5, units = "in")


tnbc_dhx9_adarp110_mb231 <-subset(immunoblot_tnbc_dhx9, Gene == "ADARp110" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_adarp110_mb231)
summary(mb231)

tnbc_dhx9_adarp110_mb231$shRNA <- factor(as.factor(tnbc_dhx9_adarp110_mb231$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb231adarp110 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_adarp110_mb231)
mb231adarp110 <- as.data.frame(mb231adarp110$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231adarp110)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231adarp110$group1 <- out$V1
mb231adarp110$group2 <- out$V2
mb231adarp110$pval <- makeStars(mb231adarp110$pval)
mb231adarp110$y.position <- max(tnbc_dhx9_adarp110_mb231$fold_change*1.1)
mb231adarp110$Cell_line = "MB231"

tnbc_dhx9_adarp110_hcc1806 <-subset(immunoblot_tnbc_dhx9, Gene == "ADARp110" & Cell_line == "HCC1806")
hcc1806 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_adarp110_hcc1806)
summary(hcc1806)

tnbc_dhx9_adarp110_hcc1806$shRNA <- factor(as.factor(tnbc_dhx9_adarp110_hcc1806$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

hcc1806adarp110 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_adarp110_hcc1806)
hcc1806adarp110 <- as.data.frame(hcc1806adarp110$shSCR[,1:4])

out <- strsplit(as.character(row.names(hcc1806adarp110)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806adarp110$group1 <- out$V1
hcc1806adarp110$group2 <- out$V2
hcc1806adarp110$pval <- makeStars(hcc1806adarp110$pval)
hcc1806adarp110$y.position <- max(tnbc_dhx9_adarp110_hcc1806$fold_change*1.1)
hcc1806adarp110$Cell_line = "HCC1806"

tnbc_dhx9_adarp110_bt549 <-subset(immunoblot_tnbc_dhx9, Gene == "ADARp110" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_adarp110_bt549)
summary(bt549)

tnbc_dhx9_adarp110_bt549$shRNA <- factor(as.factor(tnbc_dhx9_adarp110_bt549$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

bt549adarp110 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_adarp110_bt549)
bt549adarp110 <- as.data.frame(bt549adarp110$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549adarp110)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549adarp110$group1 <- out$V1
bt549adarp110$group2 <- out$V2
bt549adarp110$pval <- makeStars(bt549adarp110$pval)
bt549adarp110$y.position <- max(tnbc_dhx9_adarp110_bt549$fold_change*1.1)
bt549adarp110$Cell_line = "BT549"


tnbc_dhx9_adarp110_mb468 <-subset(immunoblot_tnbc_dhx9, Gene == "ADARp110" & Cell_line == "MB468")
mb468 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_adarp110_mb468)
summary(mb468)

tnbc_dhx9_adarp110_mb468$shRNA <- factor(as.factor(tnbc_dhx9_adarp110_mb468$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb468adarp110 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_adarp110_mb468)
mb468adarp110 <- as.data.frame(mb468adarp110$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb468adarp110)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468adarp110$group1 <- out$V1
mb468adarp110$group2 <- out$V2
mb468adarp110$pval <- makeStars(mb468adarp110$pval)
mb468adarp110$y.position <- max(tnbc_dhx9_adarp110_mb468$fold_change*1.1)
mb468adarp110$Cell_line = "MB468"

tnbcadarp110 <- rbind(mb468adarp110, mb231adarp110, hcc1806adarp110, bt549adarp110)
tnbcadarp110 <- subset(tnbcadarp110, !pval == "ns")

ggplot(subset(immunoblot_tnbc_dhx9, Gene == "ADARp110"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDHX9.3", "shDHX9.5"), labels = c("shSCR", "shDHX9-3", "shDHX9-5")) + labs(x = "", y = "Relative ADAR1-p110 Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + scale_colour_manual(values = pallete3, 
                                                                                    labels = c(BT549 = "BT-549", HCC1806  = "HCC1806", MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")) +
  geom_col(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "ADARp110", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "ADARp110", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcadarp110, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/adarp110_tnbc_dhx9_facet.tiff", height =5, width = 5, units = "in")


tnbc_dhx9_pkr_mb231 <-subset(immunoblot_tnbc_dhx9, Gene == "total PKR" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_pkr_mb231)
summary(mb231)

tnbc_dhx9_pkr_mb231$shRNA <- factor(as.factor(tnbc_dhx9_pkr_mb231$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb231pkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_pkr_mb231)
mb231pkr <- as.data.frame(mb231pkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231pkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231pkr$group1 <- out$V1
mb231pkr$group2 <- out$V2
mb231pkr$pval <- makeStars(mb231pkr$pval)
mb231pkr$y.position <- max(tnbc_dhx9_pkr_mb231$fold_change*1.1)
mb231pkr$Cell_line = "MB231"

tnbc_dhx9_pkr_hcc1806 <-subset(immunoblot_tnbc_dhx9, Gene == "total PKR" & Cell_line == "HCC1806")
hcc1806 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_pkr_hcc1806)
summary(hcc1806)

tnbc_dhx9_pkr_hcc1806$shRNA <- factor(as.factor(tnbc_dhx9_pkr_hcc1806$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

hcc1806pkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_pkr_hcc1806)
hcc1806pkr <- as.data.frame(hcc1806pkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(hcc1806pkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806pkr$group1 <- out$V1
hcc1806pkr$group2 <- out$V2
hcc1806pkr$pval <- makeStars(hcc1806pkr$pval)
hcc1806pkr$y.position <- max(tnbc_dhx9_pkr_hcc1806$fold_change*1.1)
hcc1806pkr$Cell_line = "HCC1806"

tnbc_dhx9_pkr_bt549 <-subset(immunoblot_tnbc_dhx9, Gene == "total PKR" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_pkr_bt549)
summary(bt549)

tnbc_dhx9_pkr_bt549$shRNA <- factor(as.factor(tnbc_dhx9_pkr_bt549$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

bt549pkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_pkr_bt549)
bt549pkr <- as.data.frame(bt549pkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549pkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549pkr$group1 <- out$V1
bt549pkr$group2 <- out$V2
bt549pkr$pval <- makeStars(bt549pkr$pval)
bt549pkr$y.position <- max(tnbc_dhx9_pkr_bt549$fold_change*1.1)
bt549pkr$Cell_line = "BT549"


tnbc_dhx9_pkr_mb468 <-subset(immunoblot_tnbc_dhx9, Gene == "total PKR" & Cell_line == "MB468")
mb468 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_pkr_mb468)
summary(mb468)

tnbc_dhx9_pkr_mb468$shRNA <- factor(as.factor(tnbc_dhx9_pkr_mb468$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb468pkr <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_pkr_mb468)
mb468pkr <- as.data.frame(mb468pkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb468pkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468pkr$group1 <- out$V1
mb468pkr$group2 <- out$V2
mb468pkr$pval <- makeStars(mb468pkr$pval)
mb468pkr$y.position <- max(tnbc_dhx9_pkr_mb468$fold_change*1.1)
mb468pkr$Cell_line = "MB468"

tnbcpkr <- rbind(mb468pkr, mb231pkr, hcc1806pkr, bt549pkr)
tnbcpkr <- subset(tnbcpkr, !pval == "ns")

ggplot(subset(immunoblot_tnbc_dhx9, Gene == "total PKR"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDHX9.3", "shDHX9.5"), labels = c("shSCR", "shDHX9-3", "shDHX9-5")) + labs(x = "", y = "Relative PKR Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + scale_colour_manual(values = pallete3, 
                                                                                    labels = c(BT549 = "BT-549", HCC1806  = "HCC1806", MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")) +
  geom_col(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "total PKR", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "total PKR", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcpkr, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/pkr_tnbc_dhx9_facet.tiff", height =5, width = 5, units = "in")




tnbc_dhx9_isg15_mb231 <-subset(immunoblot_tnbc_dhx9, Gene == "ISG15" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_isg15_mb231)
summary(mb231)

tnbc_dhx9_isg15_mb231$shRNA <- factor(as.factor(tnbc_dhx9_isg15_mb231$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb231isg15 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_isg15_mb231)
mb231isg15 <- as.data.frame(mb231isg15$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231isg15)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231isg15$group1 <- out$V1
mb231isg15$group2 <- out$V2
mb231isg15$pval <- makeStars(mb231isg15$pval)
mb231isg15$y.position <- max(tnbc_dhx9_isg15_mb231$fold_change*1.1)
mb231isg15$Cell_line = "MB231"

tnbc_dhx9_isg15_hcc1806 <-subset(immunoblot_tnbc_dhx9, Gene == "ISG15" & Cell_line == "HCC1806")
hcc1806 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_isg15_hcc1806)
summary(hcc1806)

tnbc_dhx9_isg15_hcc1806$shRNA <- factor(as.factor(tnbc_dhx9_isg15_hcc1806$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

hcc1806isg15 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_isg15_hcc1806)
hcc1806isg15 <- as.data.frame(hcc1806isg15$shSCR[,1:4])

out <- strsplit(as.character(row.names(hcc1806isg15)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806isg15$group1 <- out$V1
hcc1806isg15$group2 <- out$V2
hcc1806isg15$pval <- makeStars(hcc1806isg15$pval)
hcc1806isg15$y.position <- max(tnbc_dhx9_isg15_hcc1806$fold_change*1.1)
hcc1806isg15$Cell_line = "HCC1806"

tnbc_dhx9_isg15_bt549 <-subset(immunoblot_tnbc_dhx9, Gene == "ISG15" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_isg15_bt549)
summary(bt549)

tnbc_dhx9_isg15_bt549$shRNA <- factor(as.factor(tnbc_dhx9_isg15_bt549$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

bt549isg15 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_isg15_bt549)
bt549isg15 <- as.data.frame(bt549isg15$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549isg15)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549isg15$group1 <- out$V1
bt549isg15$group2 <- out$V2
bt549isg15$pval <- makeStars(bt549isg15$pval)
bt549isg15$y.position <- max(tnbc_dhx9_isg15_bt549$fold_change*1.1)
bt549isg15$Cell_line = "BT549"


tnbc_dhx9_isg15_mb468 <-subset(immunoblot_tnbc_dhx9, Gene == "ISG15" & Cell_line == "MB468")
mb468 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_isg15_mb468)
summary(mb468)

tnbc_dhx9_isg15_mb468$shRNA <- factor(as.factor(tnbc_dhx9_isg15_mb468$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb468isg15 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_isg15_mb468)
mb468isg15 <- as.data.frame(mb468isg15$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb468isg15)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468isg15$group1 <- out$V1
mb468isg15$group2 <- out$V2
mb468isg15$pval <- makeStars(mb468isg15$pval)
mb468isg15$y.position <- max(tnbc_dhx9_isg15_mb468$fold_change*1.1)
mb468isg15$Cell_line = "MB468"

tnbcisg15 <- rbind(mb468isg15, mb231isg15, hcc1806isg15, bt549isg15)
tnbcisg15 <- subset(tnbcisg15, !pval == "ns")

ggplot(subset(immunoblot_tnbc_dhx9, Gene == "ISG15"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDHX9.3", "shDHX9.5"), labels = c("shSCR", "shDHX9-3", "shDHX9-5")) + labs(x = "", y = "Relative ISG15 Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + scale_colour_manual(values = pallete3, 
                                                                                    labels = c(BT549 = "BT-549", HCC1806  = "HCC1806", MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")) +
  geom_col(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "ISG15", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "ISG15", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcisg15, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/isg15_tnbc_dhx9_facet.tiff", height =5, width = 5, units = "in")








tnbc_dhx9_dhx9_mb231 <-subset(immunoblot_tnbc_dhx9, Gene == "DHX9" & Cell_line == "MB231")
mb231 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_dhx9_mb231)
summary(mb231)

tnbc_dhx9_dhx9_mb231$shRNA <- factor(as.factor(tnbc_dhx9_dhx9_mb231$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb231dhx9 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_dhx9_mb231)
mb231dhx9 <- as.data.frame(mb231dhx9$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231dhx9)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231dhx9$group1 <- out$V1
mb231dhx9$group2 <- out$V2
mb231dhx9$pval <- makeStars(mb231dhx9$pval)
mb231dhx9$y.position <- max(tnbc_dhx9_dhx9_mb231$fold_change*1.1)
mb231dhx9$Cell_line = "MB231"

tnbc_dhx9_dhx9_hcc1806 <-subset(immunoblot_tnbc_dhx9, Gene == "DHX9" & Cell_line == "HCC1806")
hcc1806 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_dhx9_hcc1806)
summary(hcc1806)

tnbc_dhx9_dhx9_hcc1806$shRNA <- factor(as.factor(tnbc_dhx9_dhx9_hcc1806$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

hcc1806dhx9 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_dhx9_hcc1806)
hcc1806dhx9 <- as.data.frame(hcc1806dhx9$shSCR[,1:4])

out <- strsplit(as.character(row.names(hcc1806dhx9)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806dhx9$group1 <- out$V1
hcc1806dhx9$group2 <- out$V2
hcc1806dhx9$pval <- makeStars(hcc1806dhx9$pval)
hcc1806dhx9$y.position <- max(tnbc_dhx9_dhx9_hcc1806$fold_change*1.1)
hcc1806dhx9$Cell_line = "HCC1806"

tnbc_dhx9_dhx9_bt549 <-subset(immunoblot_tnbc_dhx9, Gene == "DHX9" & Cell_line == "BT549")
bt549 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_dhx9_bt549)
summary(bt549)

tnbc_dhx9_dhx9_bt549$shRNA <- factor(as.factor(tnbc_dhx9_dhx9_bt549$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

bt549dhx9 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_dhx9_bt549)
bt549dhx9 <- as.data.frame(bt549dhx9$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549dhx9)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549dhx9$group1 <- out$V1
bt549dhx9$group2 <- out$V2
bt549dhx9$pval <- makeStars(bt549dhx9$pval)
bt549dhx9$y.position <- max(tnbc_dhx9_dhx9_bt549$fold_change*1.1)
bt549dhx9$Cell_line = "BT549"


tnbc_dhx9_dhx9_mb468 <-subset(immunoblot_tnbc_dhx9, Gene == "DHX9" & Cell_line == "MB468")
mb468 <- aov(fold_change ~ shRNA, data = tnbc_dhx9_dhx9_mb468)
summary(mb468)

tnbc_dhx9_dhx9_mb468$shRNA <- factor(as.factor(tnbc_dhx9_dhx9_mb468$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb468dhx9 <- DunnettTest(fold_change ~ shRNA, data = tnbc_dhx9_dhx9_mb468)
mb468dhx9 <- as.data.frame(mb468dhx9$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb468dhx9)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468dhx9$group1 <- out$V1
mb468dhx9$group2 <- out$V2
mb468dhx9$pval <- makeStars(mb468dhx9$pval)
mb468dhx9$y.position <- max(tnbc_dhx9_dhx9_mb468$fold_change*1.1)
mb468dhx9$Cell_line = "MB468"

tnbcdhx9 <- rbind(mb468dhx9, mb231dhx9, hcc1806dhx9, bt549dhx9)
tnbcdhx9 <- subset(tnbcdhx9, !pval == "ns")

ggplot(subset(immunoblot_tnbc_dhx9, Gene == "DHX9"), aes(shRNA, fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDHX9.3", "shDHX9.5"), labels = c("shSCR", "shDHX9-3", "shDHX9-5")) + labs(x = "", y = "Relative DHX9 Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + scale_colour_manual(values = pallete3, 
                                                                                    labels = c(BT549 = "BT-549", HCC1806  = "HCC1806", MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")) +
  geom_col(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "DHX9", colour = Cell_line), aes(shRNA, mean_expression), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = subset(immunoblot_tnbc_dhx9_sum, Gene == "DHX9", colour = Cell_line), aes(x = shRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2) +
  stat_pvalue_manual(data = tnbcdhx9, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("Western_plots/dhx9_tnbc_dhx9_facet.tiff", height =5, width = 5, units = "in")


#read ff data for knockdowns
ff_bt549_dhx9 <- fread("DHX9 Foci Formation/BT549/ff_bt549_dhx9.txt")

ff_bt549_dhx9 <- na.omit(ff_bt549_dhx9)

ff_mb231_dhx9 <- fread("DHX9 Foci Formation/MB231/ff_mb231_dhx9.txt")

ff_mb231_dhx9 <- na.omit(ff_mb231_dhx9)

ff_mb468_dhx9 <- fread("DHX9 Foci Formation/MB468/ff_mb468_dhx9.txt")

ff_mb468_dhx9 <- na.omit(ff_mb468_dhx9)

ff_hcc1806_dhx9 <- fread("DHX9 Foci Formation/HCC1806/ff_hcc1806_dhx9.txt")

ff_hcc1806_dhx9 <- na.omit(ff_hcc1806_dhx9)

ff_tnbc_dhx9 <- rbind(ff_mb468_dhx9, ff_mb231_dhx9, ff_bt549_dhx9, ff_hcc1806_dhx9)

ff_tnbc_dhx9$shRNA <- gsub("-", ".", ff_tnbc_dhx9$shRNA)

#group and summarise data
ff_tnbc_dhx9_sum <- dplyr::group_by(ff_tnbc_dhx9, Cell_line, shRNA)

ff_tnbc_dhx9_sum <- dplyr::summarise(ff_tnbc_dhx9_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_tnbc_dhx9_sum$ci_area <- qnorm(0.975)*(ff_tnbc_dhx9_sum$sd_area/sqrt(3))

tnbc_dhx9_ff_mb231 <-subset(ff_tnbc_dhx9, Cell_line == "MB231")
mb231 <- aov(Relative_area ~ shRNA, data = tnbc_dhx9_ff_mb231)
summary(mb231)

tnbc_dhx9_ff_mb231$shRNA <- factor(as.factor(tnbc_dhx9_ff_mb231$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb231ff <- DunnettTest(Relative_area ~ shRNA, data = tnbc_dhx9_ff_mb231)
mb231ff <- as.data.frame(mb231ff$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb231ff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb231ff$group1 <- out$V1
mb231ff$group2 <- out$V2
mb231ff$pval <- makeStars(mb231ff$pval)
mb231ff$y.position <- max(tnbc_dhx9_ff_mb231$Relative_area*1.1)
mb231ff$Cell_line = "MB231"

tnbc_dhx9_ff_hcc1806 <-subset(ff_tnbc_dhx9, Cell_line == "HCC1806")
hcc1806 <- aov(Relative_area ~ shRNA, data = tnbc_dhx9_ff_hcc1806)
summary(hcc1806)

tnbc_dhx9_ff_hcc1806$shRNA <- factor(as.factor(tnbc_dhx9_ff_hcc1806$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

hcc1806ff <- DunnettTest(Relative_area ~ shRNA, data = tnbc_dhx9_ff_hcc1806)
hcc1806ff <- as.data.frame(hcc1806ff$shSCR[,1:4])

out <- strsplit(as.character(row.names(hcc1806ff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806ff$group1 <- out$V1
hcc1806ff$group2 <- out$V2
hcc1806ff$pval <- makeStars(hcc1806ff$pval)
hcc1806ff$y.position <- max(tnbc_dhx9_ff_hcc1806$Relative_area*1.1)
hcc1806ff$Cell_line = "HCC1806"

tnbc_dhx9_ff_bt549 <-subset(ff_tnbc_dhx9, Cell_line == "BT549")
bt549 <- aov(Relative_area ~ shRNA, data = tnbc_dhx9_ff_bt549)
summary(bt549)

tnbc_dhx9_ff_bt549$shRNA <- factor(as.factor(tnbc_dhx9_ff_bt549$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

bt549ff <- DunnettTest(Relative_area ~ shRNA, data = tnbc_dhx9_ff_bt549)
bt549ff <- as.data.frame(bt549ff$shSCR[,1:4])

out <- strsplit(as.character(row.names(bt549ff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ff$group1 <- out$V1
bt549ff$group2 <- out$V2
bt549ff$pval <- makeStars(bt549ff$pval)
bt549ff$y.position <- max(tnbc_dhx9_ff_bt549$Relative_area*1.1)
bt549ff$Cell_line = "BT549"


tnbc_dhx9_ff_mb468 <-subset(ff_tnbc_dhx9, Cell_line == "MB468")
mb468 <- aov(Relative_area ~ shRNA, data = tnbc_dhx9_ff_mb468)
summary(mb468)

tnbc_dhx9_ff_mb468$shRNA <- factor(as.factor(tnbc_dhx9_ff_mb468$shRNA), levels = c("shSCR", "shDHX9.3", "shDHX9.5"))

mb468ff <- DunnettTest(Relative_area ~ shRNA, data = tnbc_dhx9_ff_mb468)
mb468ff <- as.data.frame(mb468ff$shSCR[,1:4])

out <- strsplit(as.character(row.names(mb468ff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468ff$group1 <- out$V1
mb468ff$group2 <- out$V2
mb468ff$pval <- makeStars(mb468ff$pval)
mb468ff$y.position <- max(tnbc_dhx9_ff_mb468$Relative_area*1.1)
mb468ff$Cell_line = "MB468"

tnbcff <- rbind(mb468ff, mb231ff, hcc1806ff, bt549ff)
tnbcff <- subset(tnbcff, !pval == "ns")

ggplot(ff_tnbc_dhx9, aes(shRNA, Relative_area, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDHX9.3", "shDHX9.5"), labels = c("shSCR", "shDHX9-3", "shDHX9-5")) + labs(x = "", y = "Relative Foci Area", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(shape = "none") + facet_grid(cols = vars(Cell_line)) + scale_colour_manual(values = pallete3, 
                                                                                    labels = c(BT549 = "BT-549", HCC1806  = "HCC1806", MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")) +
  geom_col(data = ff_tnbc_dhx9_sum, aes(shRNA, mean_area, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = ff_tnbc_dhx9_sum, aes(x = shRNA, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area, colour = Cell_line), width=0.2, position = position_dodge(width = 0.7)) + 
  stat_pvalue_manual(data = tnbcff, step.increase = 0.05, step.group.by = "Cell_line", color = "Cell_line")
ggsave("FF_plots/tnbc_dhx9_facet.tiff", height =5, width = 5, units = "in")







######## This section for the rescue in SKBR3 with DHX9 or mutants
#read western (immunoblot) data for knockdowns
immunoblot_skbr3_rescue_dhx9 <- fread("immunoblot_skbr3_DHX9_rescue.txt")

immunoblot_skbr3_rescue_dhx9 <- na.omit(immunoblot_skbr3_rescue_dhx9)

immunoblot_skbr3_rescue_dhx9$sample <- paste(immunoblot_skbr3_rescue_dhx9$shRNA_2, immunoblot_skbr3_rescue_dhx9$Overexpression, sep = "/")

immunoblot_skbr3_rescue_dhx9$sample <- gsub("-", ".", immunoblot_skbr3_rescue_dhx9$sample)


#group and summarise data
immunoblot_skbr3_rescue_dhx9_sum <- dplyr::group_by(immunoblot_skbr3_rescue_dhx9, sample, Gene)

immunoblot_skbr3_rescue_dhx9_sum <- dplyr::summarise(immunoblot_skbr3_rescue_dhx9_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_skbr3_rescue_dhx9_sum$ci <- qnorm(0.975)*(immunoblot_skbr3_rescue_dhx9_sum$sd_expression/sqrt(4))

immunoblot_skbr3_rescue_dhx9$sample <- relevel(as.factor(immunoblot_skbr3_rescue_dhx9$sample), ref = "shSCR/EV")



immunoblot_skbr3_rescue_dhx9_ppkr <-subset(immunoblot_skbr3_rescue_dhx9, Gene == "p-PKR")
ppkr <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_ppkr)
summary(ppkr)

ppkr <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_ppkr)
ppkr <- as.data.frame(ppkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

ppkr$group1 <- out$V1
ppkr$group2 <- out$V2
ppkr$pval <- makeStars(ppkr$pval)
ppkr$y.position <- max(immunoblot_skbr3_rescue_dhx9_ppkr$fold_change*1.1)

ppkr <- subset(ppkr, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_dhx9, Gene == "p-PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/DHX9", "shSCR/DHX9_K417R", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/DHX9", "shDHX9.3/DHX9_K417R", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nDHX9", "shSCR\nDHX9-K417R", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nDHX9", "shDHX9-3\nDHX9-K417R", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "p-PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "p-PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = ppkr$group1, 
               xmax = ppkr$group2, 
               y.position = ppkr$y.position , label = ppkr$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_adar_dhx9_skbr3_rescue.tiff", height = 4.2, width = 6)




immunoblot_skbr3_rescue_dhx9_peif2a <-subset(immunoblot_skbr3_rescue_dhx9, Gene == "p-eIF2a")
peif2a <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_peif2a)
summary(peif2a)

peif2a <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_peif2a)
peif2a <- as.data.frame(peif2a$shSCR[,1:4])

out <- strsplit(as.character(row.names(peif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

peif2a$group1 <- out$V1
peif2a$group2 <- out$V2
peif2a$pval <- makeStars(peif2a$pval)
peif2a$y.position <- max(immunoblot_skbr3_rescue_dhx9_peif2a$fold_change*1.1)

peif2a <- subset(peif2a, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_dhx9, Gene == "p-eIF2a"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/DHX9", "shSCR/DHX9_K417R", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/DHX9", "shDHX9.3/DHX9_K417R", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nDHX9", "shSCR\nDHX9-K417R", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nDHX9", "shDHX9-3\nDHX9-K417R", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative p-eIF2a/eIF2a Abundance     ", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "p-eIF2a"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "p-eIF2a"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = peif2a$group1, 
               xmax = peif2a$group2, 
               y.position = peif2a$y.position , label = peif2a$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/peif2a_adar_dhx9_skbr3_rescue.tiff", height = 4.2, width = 6)




immunoblot_skbr3_rescue_dhx9_cparp <-subset(immunoblot_skbr3_rescue_dhx9, Gene == "c-PARP")
cparp <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_cparp)
summary(cparp)

ggplot(subset(immunoblot_skbr3_rescue_dhx9, Gene == "c-PARP"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/DHX9", "shSCR/DHX9_K417R", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/DHX9", "shDHX9.3/DHX9_K417R", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nDHX9", "shSCR\nDHX9-K417R", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nDHX9", "shDHX9-3\nDHX9-K417R", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "c-PARP"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "c-PARP"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2)
ggsave("Western_plots/cparp_adar_dhx9_skbr3_rescue.tiff", height = 4.2, width = 6)


immunoblot_skbr3_rescue_dhx9_dhx9 <-subset(immunoblot_skbr3_rescue_dhx9, Gene == "DHX9")
dhx9 <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_dhx9)
summary(dhx9)

dhx9 <- DunnettTest(x = immunoblot_skbr3_rescue_dhx9_dhx9$fold_change, g = immunoblot_skbr3_rescue_dhx9_dhx9$sample)
dhx9 <- as.data.frame(dhx9$`shSCR/EV`)

out <- strsplit(as.character(row.names(dhx9)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

dhx9$group1 <- out$V1
dhx9$group2 <- out$V2
dhx9$pval <- makeStars(dhx9$pval)
dhx9$y.position <- max(immunoblot_skbr3_rescue_dhx9_dhx9$fold_change*1.1)

dhx9 <- subset(dhx9, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_dhx9, Gene == "DHX9"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.15) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/DHX9", "shSCR/DHX9_K417R", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/DHX9", "shDHX9.3/DHX9_K417R", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nDHX9", "shSCR\nDHX9-K417R", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nDHX9", "shDHX9-3\nDHX9-K417R", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative DHX9 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "DHX9"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "DHX9"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = dhx9$group1, 
               xmax = dhx9$group2, 
               y.position = dhx9$y.position , label = dhx9$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/dhx9_adar_dhx9_skbr3_rescue.tiff", height = 4.2, width = 6)




immunoblot_skbr3_rescue_dhx9_adar <-subset(immunoblot_skbr3_rescue_dhx9, Gene == "ADARp110")
adar <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_adar)
summary(adar)

adar <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_adar)
adar <- as.data.frame(adar$`shSCR/EV`)

out <- strsplit(as.character(row.names(adar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

adar$group1 <- out$V1
adar$group2 <- out$V2
adar$pval <- makeStars(adar$pval)
adar$y.position <- max(immunoblot_skbr3_rescue_dhx9_adar$fold_change*1.1)

adar <- subset(adar, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_dhx9, Gene == "ADARp110"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/DHX9", "shSCR/DHX9_K417R", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/DHX9", "shDHX9.3/DHX9_K417R", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nDHX9", "shSCR\nDHX9-K417R", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nDHX9", "shDHX9-3\nDHX9-K417R", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative ADAR1-p110 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "ADARp110"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "ADARp110"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = adar$group1, 
               xmax = adar$group2, 
               y.position = adar$y.position , label = adar$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/adar_adar_dhx9_skbr3_rescue.tiff", height = 4.2, width = 6)



immunoblot_skbr3_rescue_dhx9_pkr <-subset(immunoblot_skbr3_rescue_dhx9, Gene == "total PKR")
pkr <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_pkr)
summary(pkr)

pkr <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_dhx9_pkr)
pkr <- as.data.frame(pkr$`shSCR/EV`)

out <- strsplit(as.character(row.names(pkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

pkr$group1 <- out$V1
pkr$group2 <- out$V2
pkr$pval <- makeStars(pkr$pval)
pkr$y.position <- max(immunoblot_skbr3_rescue_dhx9_pkr$fold_change*1.1)

pkr <- subset(pkr, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_dhx9, Gene == "total PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/DHX9", "shSCR/DHX9_K417R", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/DHX9", "shDHX9.3/DHX9_K417R", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nDHX9", "shSCR\nDHX9-K417R", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nDHX9", "shDHX9-3\nDHX9-K417R", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "total PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_dhx9_sum, Gene == "total PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = pkr$group1, 
               xmax = pkr$group2, 
               y.position = pkr$y.position , label = pkr$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/pkr_adar_dhx9_skbr3_rescue.tiff", height = 4.2, width = 6)

ppkr <- NULL
pkr <- NULL
dhx9 <- NULL
adar <- NULL
peif2a <- NULL
cparp <- NULL



#read ff data for knockdowns
ff_skbr3_rescue_dhx9 <- fread("DHX9 Foci Formation/SKBR3 Rescue/ff_skbr3_rescue_dhx9.txt")

ff_skbr3_rescue_dhx9 <- na.omit(ff_skbr3_rescue_dhx9)

ff_skbr3_rescue_dhx9$sample <- paste(ff_skbr3_rescue_dhx9$shRNA_2, ff_skbr3_rescue_dhx9$Condition, sep = "/")

ff_skbr3_rescue_dhx9$sample <- gsub("-", ".", ff_skbr3_rescue_dhx9$sample)

ff_skbr3_rescue_dhx9$sample <- relevel(as.factor(ff_skbr3_rescue_dhx9$sample), ref = "shSCR/EV")

ff_skbr3_rescue_dhx9$Relative_area <- ff_skbr3_rescue_dhx9$Relative_area*100

#group and summarise data
ff_skbr3_rescue_dhx9_sum <- dplyr::group_by(ff_skbr3_rescue_dhx9, sample)

ff_skbr3_rescue_dhx9_sum <- dplyr::summarise(ff_skbr3_rescue_dhx9_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_skbr3_rescue_dhx9_sum$ci <- qnorm(0.975)*(ff_skbr3_rescue_dhx9_sum$sd_area/sqrt(4))




aov_ff <- aov(Relative_area ~ sample, ff_skbr3_rescue_dhx9)
summary(aov_ff)

dff <- DunnettTest(Relative_area ~ sample, data = ff_skbr3_rescue_dhx9)
dff <- as.data.frame(dff$`shSCR/EV`)

out <- strsplit(as.character(row.names(dff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

dff$group1 <- out$V1
dff$group2 <- out$V2
dff$pval <- makeStars(dff$pval)
dff$y.position <- max(ff_skbr3_rescue_dhx9$Relative_area*1.1)

dff <- subset(dff, !pval == "ns")

ggplot(ff_skbr3_rescue_dhx9, aes(sample, Relative_area)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/DHX9", "shSCR/DHX9_K417R", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/DHX9", "shDHX9.3/DHX9_K417R", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nDHX9", "shSCR\nDHX9-K417R", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nDHX9", "shDHX9-3\nDHX9-K417R", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = ff_skbr3_rescue_dhx9_sum, aes(sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = ff_skbr3_rescue_dhx9_sum, aes(x = sample, y = mean_area, ymin = mean_area - sd_area, ymax = mean_area + sd_area), width=0.2) +
  geom_bracket(xmin = dff$group1, 
               xmax = dff$group2, 
               y.position = dff$y.position , label = dff$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("FF_plots/adar_dhx9_skbr3_rescue.tiff", height = 4.2, width = 6, units = "in")

aov_ff <- NULL
dff <- NULL






#####Rescue with ADAR overexpression

#read western (immunoblot) data for knockdowns
immunoblot_skbr3_rescue_adar <- fread("immunoblot_skbr3_DHX9_rescue_adar.txt")

immunoblot_skbr3_rescue_adar <- na.omit(immunoblot_skbr3_rescue_adar)

immunoblot_skbr3_rescue_adar$sample <- paste(immunoblot_skbr3_rescue_adar$shRNA_2, immunoblot_skbr3_rescue_adar$Overexpression, sep = "/")

immunoblot_skbr3_rescue_adar$sample <- gsub("-", ".", immunoblot_skbr3_rescue_adar$sample)


#group and summarise data
immunoblot_skbr3_rescue_adar_sum <- dplyr::group_by(immunoblot_skbr3_rescue_adar, sample, Gene)

immunoblot_skbr3_rescue_adar_sum <- dplyr::summarise(immunoblot_skbr3_rescue_adar_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_skbr3_rescue_adar_sum$ci <- qnorm(0.975)*(immunoblot_skbr3_rescue_adar_sum$sd_expression/sqrt(4))

immunoblot_skbr3_rescue_adar$sample <- relevel(as.factor(immunoblot_skbr3_rescue_adar$sample), ref = "shSCR/EV")



immunoblot_skbr3_rescue_adar_ppkr <-subset(immunoblot_skbr3_rescue_adar, Gene == "p-PKR")
ppkr <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_ppkr)
summary(ppkr)

ppkr <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_ppkr)
ppkr <- as.data.frame(ppkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

ppkr$group1 <- out$V1
ppkr$group2 <- out$V2
ppkr$pval <- makeStars(ppkr$pval)
ppkr$y.position <- max(immunoblot_skbr3_rescue_adar_ppkr$fold_change*1.1)

ppkr <- subset(ppkr, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_adar, Gene == "p-PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/ADAR1.p110", "shSCR/ADAR1.p150", 
                              "shDHX9.3/EV", "shDHX9.3/ADAR1.p110", "shDHX9.3/ADAR1.p150"),
                   labels = c("shSCR\nEV", "shSCR\nADAR1-p110", "shSCR\nADAR1-p150", 
                              "shDHX9-3\nEV", "shDHX9-3\nADAR1-p110", "shDHX9-3\nADAR1-p150")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "p-PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "p-PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = ppkr$group1, 
               xmax = ppkr$group2, 
               y.position = ppkr$y.position , label = ppkr$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_adar_dhx9_skbr3_rescue_adar.tiff", height = 4.2, width = 6)




immunoblot_skbr3_rescue_adar_peif2a <-subset(immunoblot_skbr3_rescue_adar, Gene == "p-eIF2a")
peif2a <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_peif2a)
summary(peif2a)

peif2a <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_peif2a)
peif2a <- as.data.frame(peif2a$shSCR[,1:4])

out <- strsplit(as.character(row.names(peif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

peif2a$group1 <- out$V1
peif2a$group2 <- out$V2
peif2a$pval <- makeStars(peif2a$pval)
peif2a$y.position <- max(immunoblot_skbr3_rescue_adar_peif2a$fold_change*1.1)

peif2a <- subset(peif2a, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_adar, Gene == "p-eIF2a"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/ADAR1.p110", "shSCR/ADAR1.p150", 
                              "shDHX9.3/EV", "shDHX9.3/ADAR1.p110", "shDHX9.3/ADAR1.p150"),
                   labels = c("shSCR\nEV", "shSCR\nADAR1-p110", "shSCR\nADAR1-p150", 
                              "shDHX9-3\nEV", "shDHX9-3\nADAR1-p110", "shDHX9-3\nADAR1-p150")) + 
  labs(x = "", y = "Relative p-eIF2a/eIF2a Abundance     ", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "p-eIF2a"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "p-eIF2a"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) 
ggsave("Western_plots/peif2a_adar_dhx9_skbr3_rescue_adar.tiff", height = 4.2, width = 6)




immunoblot_skbr3_rescue_adar_cparp <-subset(immunoblot_skbr3_rescue_adar, Gene == "c-PARP")
cparp <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_cparp)
summary(cparp)


cparp <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_cparp)
cparp <- as.data.frame(cparp$shSCR[,1:4])

out <- strsplit(as.character(row.names(cparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

cparp$group1 <- out$V1
cparp$group2 <- out$V2
cparp$pval <- makeStars(cparp$pval)
cparp$y.position <- max(immunoblot_skbr3_rescue_adar_cparp$fold_change*1.1)

cparp <- subset(cparp, !pval == "ns")



ggplot(subset(immunoblot_skbr3_rescue_adar, Gene == "c-PARP"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/ADAR1.p110", "shSCR/ADAR1.p150", 
                              "shDHX9.3/EV", "shDHX9.3/ADAR1.p110", "shDHX9.3/ADAR1.p150"),
                   labels = c("shSCR\nEV", "shSCR\nADAR1-p110", "shSCR\nADAR1-p150", 
                              "shDHX9-3\nEV", "shDHX9-3\nADAR1-p110", "shDHX9-3\nADAR1-p150")) + 
  labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "c-PARP"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "c-PARP"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = cparp$group1, 
               xmax = cparp$group2, 
               y.position = cparp$y.position , label = cparp$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/cparp_adar_dhx9_skbr3_rescue_adar.tiff", height = 4.2, width = 6)


immunoblot_skbr3_rescue_adar_dhx9 <-subset(immunoblot_skbr3_rescue_adar, Gene == "DHX9")
dhx9 <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_dhx9)
summary(dhx9)

dhx9 <- DunnettTest(x = immunoblot_skbr3_rescue_adar_dhx9$fold_change, g = immunoblot_skbr3_rescue_adar_dhx9$sample)
dhx9 <- as.data.frame(dhx9$`shSCR/EV`)

out <- strsplit(as.character(row.names(dhx9)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

dhx9$group1 <- out$V1
dhx9$group2 <- out$V2
dhx9$pval <- makeStars(dhx9$pval)
dhx9$y.position <- max(immunoblot_skbr3_rescue_adar_dhx9$fold_change*1.1)

dhx9 <- subset(dhx9, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_adar, Gene == "DHX9"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.15) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/ADAR1.p110", "shSCR/ADAR1.p150", 
                              "shDHX9.3/EV", "shDHX9.3/ADAR1.p110", "shDHX9.3/ADAR1.p150"),
                   labels = c("shSCR\nEV", "shSCR\nADAR1-p110", "shSCR\nADAR1-p150", 
                              "shDHX9-3\nEV", "shDHX9-3\nADAR1-p110", "shDHX9-3\nADAR1-p150")) + 
  labs(x = "", y = "Relative DHX9 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "DHX9"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "DHX9"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = dhx9$group1, 
               xmax = dhx9$group2, 
               y.position = dhx9$y.position , label = dhx9$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/dhx9_adar_dhx9_skbr3_rescue_adar.tiff", height = 4.2, width = 6)




immunoblot_skbr3_rescue_adar_adar <-subset(immunoblot_skbr3_rescue_adar, Gene == "ADARp110")
adar <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_adar)
summary(adar)

adar <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_adar)
adar <- as.data.frame(adar$`shSCR/EV`)

out <- strsplit(as.character(row.names(adar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

adar$group1 <- out$V1
adar$group2 <- out$V2
adar$pval <- makeStars(adar$pval)
adar$y.position <- max(immunoblot_skbr3_rescue_adar_adar$fold_change*1.1)

adar <- subset(adar, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_adar, Gene == "ADARp110"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/ADAR1.p110", "shSCR/ADAR1.p150", 
                              "shDHX9.3/EV", "shDHX9.3/ADAR1.p110", "shDHX9.3/ADAR1.p150"),
                   labels = c("shSCR\nEV", "shSCR\nADAR1-p110", "shSCR\nADAR1-p150", 
                              "shDHX9-3\nEV", "shDHX9-3\nADAR1-p110", "shDHX9-3\nADAR1-p150")) + 
  labs(x = "", y = "Relative ADAR1-p110 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "ADARp110"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "ADARp110"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = adar$group1, 
               xmax = adar$group2, 
               y.position = adar$y.position , label = adar$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/adar_adar_dhx9_skbr3_rescue_adar.tiff", height = 4.2, width = 6)

adar <- NULL

immunoblot_skbr3_rescue_adar_adarp150<-subset(immunoblot_skbr3_rescue_adar, Gene == "ADARp150")
adar <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_adarp150)
summary(adar)

adar <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_adarp150)
adar <- as.data.frame(adar$`shSCR/EV`)

out <- strsplit(as.character(row.names(adar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

adar$group1 <- out$V1
adar$group2 <- out$V2
adar$pval <- makeStars(adar$pval)
adar$y.position <- max(immunoblot_skbr3_rescue_adar_adarp150$fold_change*1.1)

adar <- subset(adar, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_adar, Gene == "ADARp150"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/ADAR1.p110", "shSCR/ADAR1.p150", 
                              "shDHX9.3/EV", "shDHX9.3/ADAR1.p110", "shDHX9.3/ADAR1.p150"),
                   labels = c("shSCR\nEV", "shSCR\nADAR1-p110", "shSCR\nADAR1-p150", 
                              "shDHX9-3\nEV", "shDHX9-3\nADAR1-p110", "shDHX9-3\nADAR1-p150")) + 
  labs(x = "", y = "Relative ADAR1-p150 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "ADARp150"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "ADARp150"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2)
ggsave("Western_plots/adarp150_adar_dhx9_skbr3_rescue_adar.tiff", height = 4.2, width = 6)



immunoblot_skbr3_rescue_adar_pkr <-subset(immunoblot_skbr3_rescue_adar, Gene == "total PKR")
pkr <- aov(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_pkr)
summary(pkr)

pkr <- DunnettTest(fold_change ~ sample, data = immunoblot_skbr3_rescue_adar_pkr)
pkr <- as.data.frame(pkr$`shSCR/EV`)

out <- strsplit(as.character(row.names(pkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

pkr$group1 <- out$V1
pkr$group2 <- out$V2
pkr$pval <- makeStars(pkr$pval)
pkr$y.position <- max(immunoblot_skbr3_rescue_adar_pkr$fold_change*1.1)

pkr <- subset(pkr, !pval == "ns")


ggplot(subset(immunoblot_skbr3_rescue_adar, Gene == "total PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/ADAR1.p110", "shSCR/ADAR1.p150", 
                              "shDHX9.3/EV", "shDHX9.3/ADAR1.p110", "shDHX9.3/ADAR1.p150"),
                   labels = c("shSCR\nEV", "shSCR\nADAR1-p110", "shSCR\nADAR1-p150", 
                              "shDHX9-3\nEV", "shDHX9-3\nADAR1-p110", "shDHX9-3\nADAR1-p150")) + 
  labs(x = "", y = "Relative PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "total PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_skbr3_rescue_adar_sum, Gene == "total PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = pkr$group1, 
               xmax = pkr$group2, 
               y.position = pkr$y.position , label = pkr$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/pkr_adar_dhx9_skbr3_rescue_adar.tiff", height = 4.2, width = 6)




#read ff data for knockdowns
ff_skbr3_rescue_adar <- fread("DHX9 Foci Formation/SKBR3 Rescue/ff_skbr3_rescue_adar.txt")

ff_skbr3_rescue_adar <- na.omit(ff_skbr3_rescue_adar)

ff_skbr3_rescue_adar$sample <- paste(ff_skbr3_rescue_adar$shRNA_2, ff_skbr3_rescue_adar$Condition, sep = "/")

ff_skbr3_rescue_adar$sample <- gsub("-", ".", ff_skbr3_rescue_adar$sample)

ff_skbr3_rescue_adar$sample <- relevel(as.factor(ff_skbr3_rescue_adar$sample), ref = "shSCR/EV")

ff_skbr3_rescue_adar$Relative_area <- ff_skbr3_rescue_adar$Relative_area*100

#group and summarise data
ff_skbr3_rescue_adar_sum <- dplyr::group_by(ff_skbr3_rescue_adar, sample)

ff_skbr3_rescue_adar_sum <- dplyr::summarise(ff_skbr3_rescue_adar_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_skbr3_rescue_adar_sum$ci <- qnorm(0.975)*(ff_skbr3_rescue_adar_sum$sd_area/sqrt(4))




aov_ff <- aov(Relative_area ~ sample, ff_skbr3_rescue_adar)
summary(aov_ff)

dff <- DunnettTest(Relative_area ~ sample, data = ff_skbr3_rescue_adar)
dff <- as.data.frame(dff$`shSCR/EV`)

out <- strsplit(as.character(row.names(dff)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

dff$group1 <- out$V1
dff$group2 <- out$V2
dff$pval <- makeStars(dff$pval)
dff$y.position <- max(ff_skbr3_rescue_adar$Relative_area*1.1)

dff <- subset(dff, !pval == "ns")

ggplot(ff_skbr3_rescue_adar, aes(sample, Relative_area)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/ADAR1.p110", "shSCR/ADAR1.p150", 
                              "shDHX9.3/EV", "shDHX9.3/ADAR1.p110", "shDHX9.3/ADAR1.p150"),
                   labels = c("shSCR\nEV", "shSCR\nADAR1-p110", "shSCR\nADAR1-p150", 
                              "shDHX9-3\nEV", "shDHX9-3\nADAR1-p110", "shDHX9-3\nADAR1-p150")) + 
  labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = ff_skbr3_rescue_adar_sum, aes(sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = ff_skbr3_rescue_adar_sum, aes(x = sample, y = mean_area, ymin = mean_area - sd_area, ymax = mean_area + sd_area), width=0.2) +
  geom_bracket(xmin = dff$group1, 
               xmax = dff$group2, 
               y.position = dff$y.position , label = dff$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("FF_plots/adar_dhx9_skbr3_rescue_adar.tiff", height = 4.2, width = 6, units = "in")





######## This section for the rescue in MCF7 with DHX9-dsRBD
#read western (immunoblot) data for knockdowns
immunoblot_mcf7_rescue_dhx9 <- fread("/kacottre/DHX9/immunoblot_mcf7_dhx9_rescue.txt")

immunoblot_mcf7_rescue_dhx9 <- na.omit(immunoblot_mcf7_rescue_dhx9)

immunoblot_mcf7_rescue_dhx9$sample <- paste(immunoblot_mcf7_rescue_dhx9$shRNA_2, immunoblot_mcf7_rescue_dhx9$Overexpression, sep = "/")

immunoblot_mcf7_rescue_dhx9$sample <- gsub("-", ".", immunoblot_mcf7_rescue_dhx9$sample)


#group and summarise data
immunoblot_mcf7_rescue_dhx9_sum <- dplyr::group_by(immunoblot_mcf7_rescue_dhx9, sample, Gene)

immunoblot_mcf7_rescue_dhx9_sum <- dplyr::summarise(immunoblot_mcf7_rescue_dhx9_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_mcf7_rescue_dhx9_sum$ci <- qnorm(0.975)*(immunoblot_mcf7_rescue_dhx9_sum$sd_expression/sqrt(4))

immunoblot_mcf7_rescue_dhx9$sample <- relevel(as.factor(immunoblot_mcf7_rescue_dhx9$sample), ref = "shSCR/EV")



immunoblot_mcf7_rescue_dhx9_ppkr <-subset(immunoblot_mcf7_rescue_dhx9, Gene == "p-PKR")
ppkr <- aov(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_ppkr)
summary(ppkr)

ppkr <- DunnettTest(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_ppkr)
ppkr <- as.data.frame(ppkr$shSCR[,1:4])

out <- strsplit(as.character(row.names(ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

ppkr$group1 <- out$V1
ppkr$group2 <- out$V2
ppkr$pval <- makeStars(ppkr$pval)
ppkr$y.position <- max(immunoblot_mcf7_rescue_dhx9_ppkr$fold_change*1.1)

ppkr <- subset(ppkr, !pval == "ns")


ggplot(subset(immunoblot_mcf7_rescue_dhx9, Gene == "p-PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "p-PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "p-PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = ppkr$group1, 
               xmax = ppkr$group2, 
               y.position = ppkr$y.position , label = ppkr$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)




immunoblot_mcf7_rescue_dhx9_peif2a <-subset(immunoblot_mcf7_rescue_dhx9, Gene == "p-eIF2a")
peif2a <- aov(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_peif2a)
summary(peif2a)

peif2a <- DunnettTest(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_peif2a)
peif2a <- as.data.frame(peif2a$shSCR[,1:4])

out <- strsplit(as.character(row.names(peif2a)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

peif2a$group1 <- out$V1
peif2a$group2 <- out$V2
peif2a$pval <- makeStars(peif2a$pval)
peif2a$y.position <- max(immunoblot_mcf7_rescue_dhx9_peif2a$fold_change*1.1)

peif2a <- subset(peif2a, !pval == "ns")


ggplot(subset(immunoblot_mcf7_rescue_dhx9, Gene == "p-eIF2a"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative p-eIF2a/eIF2a Abundance     ", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "p-eIF2a"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "p-eIF2a"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2)
ggsave("Western_plots/peif2a_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)




immunoblot_mcf7_rescue_dhx9_cparp <-subset(immunoblot_mcf7_rescue_dhx9, Gene == "c-PARP")
cparp <- aov(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_cparp)
summary(cparp)

cparp <- DunnettTest(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_cparp)
cparp <- as.data.frame(cparp$shSCR[,1:4])

out <- strsplit(as.character(row.names(cparp)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

cparp$group1 <- out$V1
cparp$group2 <- out$V2
cparp$pval <- makeStars(cparp$pval)
cparp$y.position <- max(immunoblot_mcf7_rescue_dhx9_cparp$fold_change*1.1)

cparp <- subset(cparp, !pval == "ns")




ggplot(subset(immunoblot_mcf7_rescue_dhx9, Gene == "c-PARP"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative c-PARP Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "c-PARP"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "c-PARP"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = cparp$group1, 
               xmax = cparp$group2, 
               y.position = cparp$y.position , label = cparp$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/cparp_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)


immunoblot_mcf7_rescue_dhx9_dhx9 <-subset(immunoblot_mcf7_rescue_dhx9, Gene == "DHX9")
dhx9 <- aov(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_dhx9)
summary(dhx9)

dhx9 <- DunnettTest(x = immunoblot_mcf7_rescue_dhx9_dhx9$fold_change, g = immunoblot_mcf7_rescue_dhx9_dhx9$sample)
dhx9 <- as.data.frame(dhx9$`shSCR/EV`)

out <- strsplit(as.character(row.names(dhx9)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

dhx9$group1 <- out$V1
dhx9$group2 <- out$V2
dhx9$pval <- makeStars(dhx9$pval)
dhx9$y.position <- max(immunoblot_mcf7_rescue_dhx9_dhx9$fold_change*1.1)

dhx9 <- subset(dhx9, !pval == "ns")


ggplot(subset(immunoblot_mcf7_rescue_dhx9, Gene == "DHX9"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.15) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative DHX9 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "DHX9"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "DHX9"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = dhx9$group1, 
               xmax = dhx9$group2, 
               y.position = dhx9$y.position , label = dhx9$pval, 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/dhx9_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)




immunoblot_mcf7_rescue_dhx9_adar <-subset(immunoblot_mcf7_rescue_dhx9, Gene == "ADARp110")
adar <- aov(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_adar)
summary(adar)

adar <- DunnettTest(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_adar)
adar <- as.data.frame(adar$`shSCR/EV`)

out <- strsplit(as.character(row.names(adar)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

adar$group1 <- out$V1
adar$group2 <- out$V2
adar$pval <- makeStars(adar$pval)
adar$y.position <- max(immunoblot_mcf7_rescue_dhx9_adar$fold_change*1.1)

adar <- subset(adar, !pval == "ns")


ggplot(subset(immunoblot_mcf7_rescue_dhx9, Gene == "ADARp110"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative ADAR1-p110 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "ADARp110"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "ADARp110"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2)
ggsave("Western_plots/adar_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)



immunoblot_mcf7_rescue_dhx9_pkr <-subset(immunoblot_mcf7_rescue_dhx9, Gene == "total PKR")
pkr <- aov(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_pkr)
summary(pkr)

pkr <- DunnettTest(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_pkr)
pkr <- as.data.frame(pkr$`shSCR/EV`)

out <- strsplit(as.character(row.names(pkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

pkr$group1 <- out$V1
pkr$group2 <- out$V2
pkr$pval <- makeStars(pkr$pval)
pkr$y.position <- max(immunoblot_mcf7_rescue_dhx9_pkr$fold_change*1.1)

pkr <- subset(pkr, !pval == "ns")


ggplot(subset(immunoblot_mcf7_rescue_dhx9, Gene == "total PKR"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "total PKR"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "total PKR"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2)
ggsave("Western_plots/pkr_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)



immunoblot_mcf7_rescue_dhx9_isg15 <-subset(immunoblot_mcf7_rescue_dhx9, Gene == "ISG15")
isg15 <- aov(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_isg15)
summary(isg15)

isg15 <- DunnettTest(fold_change ~ sample, data = immunoblot_mcf7_rescue_dhx9_isg15)
isg15 <- as.data.frame(isg15$shSCR[,1:4])

out <- strsplit(as.character(row.names(isg15)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

isg15$group1 <- out$V1
isg15$group2 <- out$V2
isg15$pval <- makeStars(isg15$pval)
isg15$y.position <- max(immunoblot_mcf7_rescue_dhx9_isg15$fold_change*1.1)

isg15 <- subset(isg15, !pval == "ns")


ggplot(subset(immunoblot_mcf7_rescue_dhx9, Gene == "ISG15"), aes(sample, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative ISG15 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "ISG15"), aes(sample, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mcf7_rescue_dhx9_sum, Gene == "ISG15"), aes(x = sample, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2)

ggsave("Western_plots/isg15_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)




ppkr <- NULL
pkr <- NULL
dhx9 <- NULL
adar <- NULL
peif2a <- NULL
cparp <- NULL




#read western (immunoblot) data for knockdowns
immunoblot_mb231_sgdhx9 <- fread("immunoblot_mb231_sgdhx9.txt")

immunoblot_mb231_sgdhx9 <- na.omit(immunoblot_mb231_sgdhx9)

immunoblot_mb231_sgdhx9 <- subset(immunoblot_mb231_sgdhx9, !`fold_change` == "Pending")


immunoblot_mb231_sgdhx9$fold_change <- as.numeric(immunoblot_mb231_sgdhx9$fold_change)

#group and summarise data
immunoblot_mb231_sgdhx9_sum <- dplyr::group_by(immunoblot_mb231_sgdhx9, sgRNA, Gene)

immunoblot_mb231_sgdhx9_sum <- dplyr::summarise(immunoblot_mb231_sgdhx9_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

immunoblot_mb231_sgdhx9_sum$ci <- qnorm(0.975)*(immunoblot_mb231_sgdhx9_sum$sd_expression/sqrt(4))

#make plots and save

# aov_ppkr <- aov(fold_change ~ sample, subset(immunoblot_mb231_sgdhx9, Gene == "p-PKR"))
# tppkr <- TukeyHSD(aov_ppkr)
# tppkr <- as.data.frame(tppkr$sample[,1:4])

immunoblot_mb231_sgdhx9$sgRNA <- factor(as.factor(immunoblot_mb231_sgdhx9$sgRNA), levels = c("NT", "sgDHX9"))

immunoblot_mb231_sgdhx9_ppkr <- subset(immunoblot_mb231_sgdhx9, Gene == "p-PKR")

ppkr <- DunnettTest(fold_change ~ sgRNA, data = immunoblot_mb231_sgdhx9_ppkr)

ggplot(subset(immunoblot_mb231_sgdhx9, Gene == "p-PKR"), aes(sgRNA, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("NT", "sgDHX9")) + 
  labs(x = "", y = "Relative p-PKR/PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mb231_sgdhx9_sum, Gene == "p-PKR"), aes(sgRNA, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mb231_sgdhx9_sum, Gene == "p-PKR"), aes(x = sgRNA, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2) +
  geom_bracket(xmin = "NT", 
               xmax = "sgDHX9", 
               y.position = max(immunoblot_mb231_sgdhx9$fold_change*1.1) , label = makeStars(ppkr$NT[4]), 
               tip.length = 0.05, step.increase = 0.1)
ggsave("Western_plots/ppkr_mb231_sgDHX9.tiff", height = 4.2, width = 4.2)



ggplot(subset(immunoblot_mb231_sgdhx9, Gene == "DHX9"), aes(sgRNA, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("NT", "sgDHX9")) + 
  labs(x = "", y = "Relative DHX9 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mb231_sgdhx9_sum, Gene == "DHX9"), aes(sgRNA, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mb231_sgdhx9_sum, Gene == "DHX9"), aes(x = sgRNA, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2)
ggsave("Western_plots/DHX9_mb231_sgDHX9.tiff", height = 4.2, width = 4.2)


ggplot(subset(immunoblot_mb231_sgdhx9, Gene == "PKR"), aes(sgRNA, fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("NT", "sgDHX9")) + 
  labs(x = "", y = "Relative PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = subset(immunoblot_mb231_sgdhx9_sum, Gene == "PKR"), aes(sgRNA, mean_expression), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = subset(immunoblot_mb231_sgdhx9_sum, Gene == "PKR"), aes(x = sgRNA, y = mean_expression, ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), width=0.2)
ggsave("Western_plots/pkr_mb231_sgDHX9.tiff", height = 4.2, width = 4.2)



######## This section for the rescue in MCF7 with DHX9-dsRBD
#read qpcr data for knockdowns
qpcr_mcf7_rescue_dhx9 <- fread("/kacottre/Weber_Lab/APEX-ADAR_manuscript_dhx9/Source Data/qPCR/qPCR_summary.txt")

qpcr_mcf7_rescue_dhx9 <- na.omit(qpcr_mcf7_rescue_dhx9)

qpcr_mcf7_rescue_dhx9$sample <- paste(qpcr_mcf7_rescue_dhx9$shRNA_2, qpcr_mcf7_rescue_dhx9$Overexpression, sep = "/")

qpcr_mcf7_rescue_dhx9$sample <- gsub("-", ".", qpcr_mcf7_rescue_dhx9$sample)


#group and summarise data
qpcr_mcf7_rescue_dhx9_sum <- dplyr::group_by(qpcr_mcf7_rescue_dhx9, sample, Pathway, Gene)

qpcr_mcf7_rescue_dhx9_sum <- dplyr::summarise(qpcr_mcf7_rescue_dhx9_sum, mean_expression = mean(fold_change), sd_expression = sd(fold_change))

qpcr_mcf7_rescue_dhx9_sum$ci <- qnorm(0.975)*(qpcr_mcf7_rescue_dhx9_sum$sd_expression/sqrt(4))



qpcr_mcf7_rescue_dhx9$sample <- relevel(as.factor(qpcr_mcf7_rescue_dhx9$sample), ref = "shSCR/GFP")

aov_qpcr <- aov(fold_change ~ sample, data = qpcr_mcf7_rescue_dhx9)
summary(aov_qpcr)

qpcr_d <- DunnettTest(fold_change ~ sample, data = qpcr_mcf7_rescue_dhx9)
qpcr_d <- as.data.frame(qpcr_d$shSCR[,1:4])

out <- strsplit(as.character(row.names(qpcr_d)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

qpcr_d$group1 <- out$V1
qpcr_d$group2 <- out$V2
qpcr_d$pval <- makeStars(qpcr_d$pval)
qpcr_d$y.position <- max(qpcr_mcf7_rescue_dhx9$fold_change*1.1)

qpcr_d <- subset(qpcr_d, !pval == "ns")

qpcr_d$Pathway = NA
qpcr_d$Gene = NA


ggplot(qpcr_mcf7_rescue_dhx9, aes(sample, fold_change, shape = Pathway, colour = Gene)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative mRNA Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette[c(1:4,6,7)])+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0, colour = "grey", linetype = "dashed") +
  guides() + stat_summary(fun = median, fun.min = median, fun.max = median,
                          geom = "crossbar", position = position_dodge(), width = 0.7) +
  geom_bracket(xmin = qpcr_d$group1, 
             xmax = qpcr_d$group2, 
             y.position = qpcr_d$y.position , label = qpcr_d$pval, 
             tip.length = 0.05, step.increase = 0.1, inherit.aes = FALSE)
ggsave("qPCR_plots/all_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6.6)


qpcr_mcf7_rescue_dhx9_isg <- subset(qpcr_mcf7_rescue_dhx9, Pathway == "ISG")

aov_qpcr <- aov(fold_change ~ sample, data = qpcr_mcf7_rescue_dhx9_isg)
summary(aov_qpcr)

qpcr_d <- DunnettTest(fold_change ~ sample, data = qpcr_mcf7_rescue_dhx9_isg)
qpcr_d <- as.data.frame(qpcr_d$shSCR[,1:4])

out <- strsplit(as.character(row.names(qpcr_d)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

qpcr_d$group1 <- out$V1
qpcr_d$group2 <- out$V2
qpcr_d$pval <- makeStars(qpcr_d$pval)
qpcr_d$y.position <- max(qpcr_mcf7_rescue_dhx9_isg$fold_change*1.1)

qpcr_d <- subset(qpcr_d, !pval == "ns")

qpcr_d$Pathway = NA
qpcr_d$Gene = NA
qpcr_d$Replicate = NA


ggplot(qpcr_mcf7_rescue_dhx9_isg, aes(sample, fold_change, colour = Gene, group = Gene)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), aes(shape = as.factor(Replicate))) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/GFP", "shDHX9.3/dsRBm.GFP"),
                   labels = c("shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative mRNA Expression", colour = "", fill = "", title = "ISGs") + 
  scale_colour_manual(values = palette[c(2,7)]) +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + stat_summary(fun = median, fun.min = median, fun.max = median,
                          geom = "crossbar", position = position_dodge(), width = 0.7) +
  geom_bracket(xmin = qpcr_d$group1, 
               xmax = qpcr_d$group2, 
               y.position = qpcr_d$y.position , label = qpcr_d$pval, 
               tip.length = 0.05, step.increase = 0.1, inherit.aes = FALSE)
ggsave("qPCR_plots/isg_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)

qpcr_d <- NULL
aov_qpcr <- NULL

qpcr_mcf7_rescue_dhx9_ATF4 <- subset(qpcr_mcf7_rescue_dhx9, Pathway == "ATF4")

aov_qpcr <- aov(fold_change ~ sample, data = qpcr_mcf7_rescue_dhx9_ATF4)
summary(aov_qpcr)

qpcr_d <- DunnettTest(fold_change ~ sample, data = qpcr_mcf7_rescue_dhx9_ATF4)
qpcr_d <- as.data.frame(qpcr_d$shSCR[,1:4])

out <- strsplit(as.character(row.names(qpcr_d)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

qpcr_d$group1 <- out$V1
qpcr_d$group2 <- out$V2
qpcr_d$pval <- makeStars(qpcr_d$pval)
qpcr_d$y.position <- max(qpcr_mcf7_rescue_dhx9_ATF4$fold_change*1.1)

qpcr_d <- subset(qpcr_d, !pval == "ns")

qpcr_d$Pathway = NA
qpcr_d$Gene = NA


ggplot(qpcr_mcf7_rescue_dhx9_ATF4, aes(sample, fold_change, colour = Gene, group = Gene)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), aes(shape = as.factor(Replicate))) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative mRNA Expression", colour = "", fill = "", title = "ATF4 Targets") + 
  scale_colour_manual(values = palette[c(1,5)]) +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + stat_summary(fun = median, fun.min = median, fun.max = median,
                          geom = "crossbar", position = position_dodge(), width = 0.7) +
  geom_bracket(xmin = qpcr_d$group1, 
               xmax = qpcr_d$group2, 
               y.position = qpcr_d$y.position , label = qpcr_d$pval, 
               tip.length = 0.05, step.increase = 0.1, inherit.aes = FALSE)
ggsave("qPCR_plots/atf4_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)

qpcr_d <- NULL
aov_qpcr <- NULL

qpcr_mcf7_rescue_dhx9_NFKB <- subset(qpcr_mcf7_rescue_dhx9, Pathway == "NF-KB")

aov_qpcr <- aov(fold_change ~ sample, data = qpcr_mcf7_rescue_dhx9_NFKB)
summary(aov_qpcr)

qpcr_d <- DunnettTest(fold_change ~ sample, data = qpcr_mcf7_rescue_dhx9_NFKB)
qpcr_d <- as.data.frame(qpcr_d$shSCR[,1:4])

out <- strsplit(as.character(row.names(qpcr_d)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

qpcr_d$group1 <- out$V1
qpcr_d$group2 <- out$V2
qpcr_d$pval <- makeStars(qpcr_d$pval)
qpcr_d$y.position <- max(qpcr_mcf7_rescue_dhx9_NFKB$fold_change*1.1)

qpcr_d <- subset(qpcr_d, !pval == "ns")

qpcr_d$Pathway = NA
qpcr_d$Gene = NA


ggplot(qpcr_mcf7_rescue_dhx9_NFKB, aes(sample, fold_change, colour = Gene, group = Gene)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), aes(shape = as.factor(Replicate))) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative mRNA Expression", colour = "", fill = "", title = "NF-KB Targets") + 
  scale_colour_manual(values = palette[c(3,6)]) +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + stat_summary(fun = median, fun.min = median, fun.max = median,
                          geom = "crossbar", position = position_dodge(), width = 0.7) +
  geom_bracket(xmin = qpcr_d$group1, 
               xmax = qpcr_d$group2, 
               y.position = qpcr_d$y.position , label = qpcr_d$pval, 
               tip.length = 0.05, step.increase = 0.1, inherit.aes = FALSE)
ggsave("qPCR_plots/nfkb_adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6)




###IF DATA
replicate_1 <- fread("replicate_1_IF.csv")
replicate_2 <- fread("replicate_2_IF.csv")
replicate_3 <- fread("replicate_3_IF.csv")
replicate_4 <- fread("replicate_4_IF.csv")

replicate_1$Replicate = 1 
replicate_2$Replicate = 2 
replicate_3$Replicate = 3 
replicate_4$Replicate = 4 

IF_data <- rbind(replicate_1, replicate_2, replicate_3, replicate_4)

IF_data$Nuclear_normalized <- IF_data$Median - IF_data$Background

IF_data$sample <- paste(IF_data$shRNA_1, IF_data$shRNA_2, sep = "/")

IF_data$sample <- gsub("-", ".", IF_data$sample)


#group and summarise data
IF_data_sum <- dplyr::group_by(IF_data, sample, Replicate)

IF_data_sum <- dplyr::summarise(IF_data_sum, mean_intensity = mean(Nuclear_normalized))

IF_data_sum_sum <- dplyr::group_by(IF_data_sum, sample)

IF_data_sum_sum <- dplyr::summarise(IF_data_sum_sum, mean_mean_intensity = mean(mean_intensity), sd_mean_intensity = sd(mean_intensity))



ggplot(IF_data, aes(sample, Nuclear_normalized, colour = as.factor(Replicate))) + geom_beeswarm(alpha = 0.1) + 
  geom_point(data = IF_data_sum, aes(sample, mean_intensity, colour = as.factor(Replicate))) +
  scale_x_discrete(limits = c("shADAR/shSCR", "shADAR/shDHX9"), 
                   labels = c("shADAR\nshSCR", "shADAR\nshDHX9")) + 
  scale_colour_manual(values = palette[1:4]) + labs(y = "Normalized Nulear J2 Intensity") +
  stat_summary(data = IF_data_sum, aes(sample, mean_intensity), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", colour = "black") + theme_science() + 
  theme(axis.title.x=element_blank(), legend.position = "none")
ggsave("J2_IF.tiff", units = "in", heigh = 4, width  = 4)

t.test(IF_data_sum$mean_intensity[IF_data_sum$sample == "shADAR/shSCR"], IF_data_sum$mean_intensity[IF_data_sum$sample == "shADAR/shDHX9"])



#read ff data for knockdowns
ff_mcf7_rescue_dhx9 <- fread("DHX9 Foci Formation/ff_mcf7_rescue_dhx9.txt")

ff_mcf7_rescue_dhx9 <- na.omit(ff_mcf7_rescue_dhx9)

ff_mcf7_rescue_dhx9$sample <- paste(ff_mcf7_rescue_dhx9$shRNA_2, ff_mcf7_rescue_dhx9$Condition, sep = "/")

ff_mcf7_rescue_dhx9$sample <- gsub("-", ".", ff_mcf7_rescue_dhx9$sample)

ff_mcf7_rescue_dhx9$sample <- relevel(as.factor(ff_mcf7_rescue_dhx9$sample), ref = "shSCR/EV")

ff_mcf7_rescue_dhx9$Relative_area <- ff_mcf7_rescue_dhx9$Relative_area *100

#group and summarise data
ff_mcf7_rescue_dhx9_sum <- dplyr::group_by(ff_mcf7_rescue_dhx9, sample)

ff_mcf7_rescue_dhx9_sum <- dplyr::summarise(ff_mcf7_rescue_dhx9_sum, mean_area = mean(Relative_area), sd_area = sd(Relative_area))

ff_mcf7_rescue_dhx9_sum$ci <- qnorm(0.975)*(ff_mcf7_rescue_dhx9_sum$sd_area/sqrt(4))



ggplot(ff_mcf7_rescue_dhx9, aes(sample, Relative_area)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR/EV", "shSCR/GFP", "shSCR/dsRBM.GFP", 
                              "shDHX9.3/EV", "shDHX9.3/GFP", "shDHX9.3/dsRBM.GFP"),
                   labels = c("shSCR\nEV", "shSCR\nGFP", "shSCR\ndsRBD-GFP", 
                              "shDHX9-3\nEV", "shDHX9-3\nGFP", "shDHX9-3\ndsRBD-GFP")) + 
  labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 100, size = 1, colour = "grey", linetype = "dashed") +
  guides(shape = "none") + 
  geom_col(data = ff_mcf7_rescue_dhx9_sum, aes(sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") +
  geom_errorbar(data = ff_mcf7_rescue_dhx9_sum, aes(x = sample, y = mean_area, ymin = mean_area - sd_area, ymax = mean_area + sd_area), width=0.2)
ggsave("FF_plots/adar_dhx9_mcf7_rescue.tiff", height = 4.2, width = 6, units = "in")





#read cell_proliferation data for knockdowns
cell_proliferation_mb231_dhx9 <- fread("Cell_proliferation_revision.txt")

cell_proliferation_mb231_dhx9 <- na.omit(cell_proliferation_mb231_dhx9)

#group and summarise data
cell_proliferation_mb231_dhx9_sum <- dplyr::group_by(cell_proliferation_mb231_dhx9, shRNA, Day)

cell_proliferation_mb231_dhx9_sum <- dplyr::summarise(cell_proliferation_mb231_dhx9_sum, mean_count = mean(Count), mean_viability = mean(Viability), 
                                                             sd_count = sd(Count), sd_viability = sd(Viability))

ggplot(cell_proliferation_mb231_dhx9_sum, aes(Day, mean_count, colour = shRNA)) +
  geom_line() + theme_science(base_size = 13) + scale_color_manual(values = palette2) +
  labs(x = "Days after transduction", y = "Viable Cell Count", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(shape = "none") + 
  geom_errorbar(aes(x = Day, y = mean_count, ymin = mean_count - sd_count, ymax = mean_count + sd_count), width=0.2) +
  geom_point(data =  cell_proliferation_mb231_dhx9, aes(Day, Count, colour = shRNA, shape = as.factor(Replicate)), alpha = 0.5)
ggsave("cell_proliferation_mb231.tiff", height = 3, width = 4.5, units = "in")



ggplot(cell_proliferation_mb231_dhx9_sum, aes(Day, mean_viability, colour = shRNA)) +
  geom_line() + theme_science(base_size = 13) + scale_color_manual(values = palette2) +
  labs(x = "Days after transduction", y = "Percent Viability", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(shape = "none") + scale_y_continuous(limits = c(0,105), breaks = c(0,20,40,60,80,100)) +
  geom_errorbar(aes(x = Day, y = mean_viability, ymin = mean_viability - sd_viability, ymax = mean_viability + sd_viability), width=0.2) +
  geom_point(data =  cell_proliferation_mb231_dhx9, aes(Day, Viability, colour = shRNA, shape = as.factor(Replicate)), alpha = 0.5)
ggsave("cell_viability_mb231.tiff", height = 3, width = 4.5, units = "in")




#read Real-Time-Glo data for knockdowns
realtime_glo_mb231_dhx9 <- fread("RealTime-Glo_revision.txt")

realtime_glo_mb231_dhx9 <- na.omit(realtime_glo_mb231_dhx9)

#group and summarise data
realtime_glo_mb231_dhx9_sum <- dplyr::group_by(realtime_glo_mb231_dhx9, shRNA, Assay)

realtime_glo_mb231_dhx9_sum <- dplyr::summarise(realtime_glo_mb231_dhx9_sum, mean_RFU_RLU = mean(RFU_RLU), sd_RFU_RLU = sd(RFU_RLU))

ggplot(realtime_glo_mb231_dhx9, aes(shRNA, RFU_RLU, colour = Assay)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shDHX9"),
                   labels = c("shSCR", "shDHX9")) + 
  labs(x = "", y = "RLU (Annexin V) or\n RFU (Membrane Integrity)", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.position = "bottom") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), trans = "log10") +
  guides(shape = "none") + scale_color_manual(values = palette2) +
  geom_col(data = realtime_glo_mb231_dhx9_sum, aes(shRNA, mean_RFU_RLU, colour = Assay), position = position_dodge(width = 0.7), width = 0.5, fill = "white", alpha = 0) +
  geom_errorbar(data = realtime_glo_mb231_dhx9_sum, aes(x = shRNA, y = mean_RFU_RLU, ymin = mean_RFU_RLU - sd_RFU_RLU, ymax = mean_RFU_RLU + sd_RFU_RLU, colour = Assay), width=0.2, position  = position_dodge(width = 0.7))
ggsave("realtime_glo_mb231.tiff", height = 4.2, width = 6, units = "in")


