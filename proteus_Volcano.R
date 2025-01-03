#read in libraries
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)
library(ggrepel)
library(svglite)
library(ggrepel)
library(tibble)

setwd("C:/Users/Olivia.Schwartz/OneDrive - University of Denver/Projects/Proteus Proteobactin/mzMINE/pos")

# Files
dset<-read.csv("Normalised_Quant_table_pos.csv", header = TRUE)
metadata <- read.csv("proteus_pos_md.csv", header = T, check.names = F, sep = ",")
colnames(metadata)[1] <- "filename"
colnames(metadata)<-gsub("ATTRIBUTE_","",colnames(metadata)) #Remove "ATTRIBUTE_"

####file names start with numbers and are therefore having an X added to the beginning upon being converted to a matrix 

## Inputs
# isotope = "18O"
# week = "8"
# strain = "WT"
colname_feat = "producer"
neg_fc_label = "no" #tg
pos_fc_label = "yes" #wt
volcano_title = "producer"

data<- as.matrix(dset[-c(1:3)])
featureIDs <- dset[,1]
# Check the loaded dataset
# Dimension of the dataset
dim(data) 
head(data)

# Separate the two conditions into two smaller data frames

filtered_metadata <- metadata %>% dplyr::filter(metadata[["bio_rep"]] == "3")


wtMetadata <- filtered_metadata[which(filtered_metadata[[colname_feat]]== paste0(pos_fc_label)), ]

wtFileNames <- c(wtMetadata$filename)
wtDSet <- dset[, wtFileNames, drop = FALSE]
wtData<- as.matrix(wtDSet %>% select(-1))

tgMetadata <- filtered_metadata[ which(filtered_metadata[[colname_feat]]== paste0(neg_fc_label)), ]

tgFileNames <- c(tgMetadata$filename)
tgDSet <- dset[, tgFileNames, drop = FALSE]
tgData<- as.matrix(tgDSet %>% select(-1))

data <- as.matrix(cbind(wtDSet,tgDSet))
featureIDs <- dset[,1]


# Check the loaded dataset
# Dimension of the dataset
dim(data) 
head(data)

# Separate the two conditions into two smaller data frames
wt = wtData
ko = tgData

# Compute the means of the samples of each condition
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)

# Just get the maximum of all the means for plotting
limit = max(wt.mean, ko.mean)
means <- cbind.data.frame(wt.mean, ko.mean)

###############################################################
######### WT VS KO ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = wt.mean / ko.mean
log2fold <- log2(fold)
#add to results dataframe
DEres.df$logFC <- log2fold

# Generate histogram of the fold differences
foldhist <- ggplot(DEres.df, aes(x=logFC)) +
  geom_histogram( binwidth=0.01, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Fold Change Distribution") +
  labs(x= "Log fold change in abundance between groups", y = "Count") +
  theme_minimal() 
foldhist

#compute stat sig - t-test
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics
for(i in 1 : nrow(data)) { # For each gene : 
  x = wt[i,] # WT of gene number i
  y = ko[i,] # KO of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y, na.rm=TRUE, paired=FALSE, exact=FALSE)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}
# add results to results dataframe
DEres.df$pvalue <- pvalue

# Generate -log10(pvalues)
pvallogged <- -log10(pvalue)
# add results to results dataframe
DEres.df$logpval <- pvallogged
# Generate histogram of -log10(pvalues)
pval_hist <- ggplot(DEres.df, aes(x=logpval)) +
  geom_histogram( binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("-log10(pvalue) Distribution") +
  labs(x= "Significance levels occuring in D.E. results", y = "Count") +
  xlim(0,4) +
  theme_minimal() 
pval_hist

DEres.df <- DEres.df %>%
  separate_wider_delim(featureIDs, delim = "_", names = c("row_ID", "row_mz","row_rt"),cols_remove = FALSE)

# Volcano Plots:
# put the biological significance (fold changes)
# and statistical significance (p-value) in one plot
# Generate the volcano plot
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          paste0(volcano_title)){
  # Inputs: ##############################################
  
  # results.df: dataframe containing columns: 
  #             1) "featureIDs": feature IDs
  #             2) "logFC": foldchange values from DE analysis
  #             3) "logpval": log-transformed p-values from DE analysis
  # fold cutoff: significance threshold to render visually on plot;
  #               denotes fold difference between mutant and wildtype
  #               also referred to as "biologial signal"
  # p_value_cutoff: significance threshold to render visually on plot;
  #                 denotes statistical significance
  ########################################################
  # create factor denoting differential expression status
  stats.df <- results.df %>% mutate(diffexpressed = 
                                      ifelse(logFC < -(fold_cutoff) & 
                                               logpval >= -log10(pvalue_cutoff),
                                             paste0(neg_fc_label), 
                                             ifelse(logFC > fold_cutoff & 
                                                      logpval >= -log10(pvalue_cutoff), 
                                                    paste0(pos_fc_label), "Not significant")))
  
  stats.df$diffexpressed <- as.factor(stats.df$diffexpressed)
  # Base plot
  plot <- ggplot(stats.df, aes(x = logFC, y = logpval, col = diffexpressed)) +
    geom_point(aes(alpha=0.5)) +
    geom_text(data=subset(stats.df, logFC < -(fold_cutoff) & 
                            logpval >= -log10(pvalue_cutoff) | logFC > fold_cutoff & 
                            logpval >= -log10(pvalue_cutoff)),
              aes(logFC,logpval,label=row_ID)) +
    # geom_label_repel() +
    ylab("-log10(pvalue)") +
    geom_vline(xintercept=c(-fold_cutoff, fold_cutoff), col="red") +
    geom_hline(yintercept=-log10(pvalue_cutoff), col="darkturquoise") +
    scale_color_manual(values=c("black", "blue", "red"), 
                       name = "Differential expression\nstatus") +
    theme_minimal()
  return(plot)
  
}

library("plotly")
# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01) + ggtitle(volcano_title)
ggplotly(plot)
# ggsave(file=paste0(volcano_title,".svg"), plot=last_plot())

# install.packages('svglite')
ggsave(file=paste0(volcano_title,".svg"), plot=plot, width=10, height=8)


sig_values <- subset(DEres.df, logFC < -0.5 & 
              logpval >= -log10(0.01) | logFC > 0.5 & 
              logpval >= -log10(0.01))

write.csv(sig_values, file=paste0(volcano_title, ".csv"))
# write.csv(DEres.df, "~/Desktop/volcano_t-test_pos.csv", row.names=FALSE)


data_anova <- t(DEres.df)
data_anova <- as.data.frame(data_anova)
colnames(data_anova) <- data_anova[4,]
data_anova <- data_anova[-(35:37),]
data_anova <- data_anova[-(1:4),]  
data_anova <- add_column(data_anova, filename = row.names(data_anova), .before = 1)

data_merge <- filtered_metadata %>%
  inner_join(.,data_anova)

# 
# stats.df, logpval >= -log10(pvalue_cutoff) & logFC > fold_cutoff | logpval >= 2.5 & logFC < fold_cutoff)
# 
# res.aov <- aov(featureIDs ~ strain_treatment, data = my_data)
# 


