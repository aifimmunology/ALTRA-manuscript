##############################################################################
#
# Script: IC_flow_bcell_data_analysis_final.R
# Project: He, Glass et al pre-clinical RA study 
# Subproject: Peripheral B cell intracellular flow cytometry data analysis
# Author: Marla Glass
# Date: 09-04-2025
#
# This program takes in a CSV file with pre-processed, isotype-labeled, and metaclustered
# B cell intracellular flow cytometry data for ACPA- HC2 (control) and ACPA+ ARI subjects (bcell_clustered_flow_data.csv)
#
# Analyses:
# Cytokine-positive B cells by status/group - Main Figure 3I analysis plot
#
#######################################################################################################


###### LIBRARIES ######

require(tidyverse)
require(ggplot2)
require(data.table)
require(viridis)
require(PNWColors)
require(reshape2)
require(ggrepel)
require(ggsignif)
require(magrittr)

###### USER INPUTS ######

# *update these file paths as needed*
images.path <- "~/data_analysis/images/"
path <- "~/data_analysis/tables/"
dat.path <- "~/data_analysis/tables/bcell_clustered_flow_data.csv"

factors <- c("cell_type", "condition", "donor", "time", "sample_status", 
             "isotype", "meta", "status", "batch")

conds <- c("CpG_CD40_stimulated", "unstimulated")
#cond.colors <- c('#2F4858', '#EFBFBD') %>% setNames(conds)

donors <- c("ARI_1" , "ARI_2" , "ARI_3" , "ARI_4" , "ARI_5" , "ARI_6" , 
            "ARI_7" , "ARI_8" , "ARI_9" , "ARI_10" , "ARI_11" , "ARI_12" , 
            "ARI_13" , "ARI_14" , "ARI_15" , "ARI_16" , "ARI_17" , 
            "HC_1" , "HC_2" , "HC_3" , "HC_4" , "HC_5" , "HC_6" , 
            "HC_7" , "HC_8" , "HC_9" , "HC_10" , "HC_11" , "HC_12" , 
            "HC_13")

batches <- c("1", "2", "3")

clin.status <- c("ARI", "HC2")
clin.colors <- c("#F59F00", "#5AAA46") %>% setNames(clin.status)

cytokines <- c("TNFa", "IL_6", "IL_10")
cytok.factors <- c("IL_10_pos", "IL_6_pos", "TNFa_pos")
#cytok.colors <- pnw_palette("Bay", 3) %>% as.vector() %>% setNames(cytok.factors)

intracellular <- c("TNFa", "IL_6", "IL_10", "RANKL")
intra.factors <- c("IL_10_pos", "IL_6_pos", "TNFa_pos", "RANKL_pos")
#intra.colors <- pnw_palette("Bay", 4) %>% as.vector() %>% setNames(intra.factors)

all.cols <- c("CD107a",	"CD69",	"CD45RB",	
              "CD24", "CD20",	"HLA_DR",	"CD27",	"CD10",	
              "CD11c", "CD73", "CD95", "CD11c", "CD21",
              "IgG", "IgA", "IgD", "IgM", 
              "TNFa", "IL_6", "IL_10", "RANKL")

subsets <- c("CD27_neg_Effector", "CD27_pos_Effector", "Early_Memory", "Core_Memory", "Late_Memory", "CD95_Memory", 
             "Transitional", "Naive", "Plasma")
#subset.colors <- c("#F9B5AC", "#861D31", "#C4E2E1", "#399390", "#849324", "#263FA6", "#805D93", "#CF7C63", "#323949") %>% setNames(subsets)

isotypes <- c("IgD", "IgMD", "IgM", "IgG", "IgA", "ND", "surface_Ig-")
#isotype.colors <- c("#664A5B", "#C05746", "#4D7184", "#8DB979", "#012A36", "#D1D2D4", "#ED7D3B") %>% setNames(isotypes)

###### FUNCTIONS ######

naiveBCytokineStatus <- function(dt, 
                            sc=clin.colors,
                            pa=images.path) {
  # Generates plots with cytokine-positive naive B cell percentages by status/groups 
  # Applies Wilcoxon rank sum test to determine statistical significance between groups
  # Inputs:
  #   dt - data.table
  #   sc - named vector of status colors
  #   pa - path to images folder
  # Outputs:
  #  Plots of cytokine-positive naive B cell %s by status/group with adjusted p values reported
  
  #TNFa+ naive B cells
  perTNFastat <- table(dt[meta=="Naive", .(donor, TNFa_pos, status)]) %>%
    as.data.table() %>%
    .[TNFa_pos==T, NTNFa:=N, by=donor] %>%
    .[, PTNFa:=100*NTNFa/sum(N), by=donor] %>%
    .[, c("N", "NTNFa"):=NULL] %>%
    .[TNFa_pos==T] %>%
    .[, TNFa_pos:=NULL] %>% 
    .[PTNFa!=0.000000] %>%
    .[, donor:=factor(donor, levels=donors)] %>%
    .[, status:=factor(status, levels=clin.status)]
  
  gt.vec <- perTNFastat[status=="ARI", PTNFa]
  tc.vec <- perTNFastat[status=="HC2", PTNFa]
  stat.test <- wilcox.test(gt.vec, tc.vec, paired=FALSE, alternative="g")
  
  ggplot(perTNFastat, aes(status, PTNFa, fill=status))  + 
    geom_boxplot(width=0.5, color="black", position="dodge") + 
    geom_jitter(aes(status, PTNFa), size=8, height=0, width=0.03) +
    ylab("TNFa+ Naive B cells (%)") +
    ggtitle("TNFa+ Naive B cells") + 
    scale_fill_manual(values=sc) +
    theme_bw() +
    theme(text=element_text(size=42), 
          title=element_text(size=40),
          legend.position='none')
  ggsave(paste0(pa, "Perc_TNFa_naiveb_status_stats.png"), width=8, height=12)
  
  #IL-6+ naive B cells
  perIL_6stat <- table(dt[meta=="Naive", .(donor, IL_6_pos, status)]) %>%
    as.data.table() %>%
    .[IL_6_pos==T, NIL_6:=N, by=donor] %>%
    .[, PIL_6:=100*NIL_6/sum(N), by=donor] %>%
    .[, c("N", "NIL_6"):=NULL] %>%
    .[IL_6_pos==T] %>%
    .[, IL_6_pos:=NULL] %>% 
    .[PIL_6!=0.000000] %>%
    .[, donor:=factor(donor, levels=donors)] %>%
    .[, status:=factor(status, levels=clin.status)]
  
  gt.vec <- perIL_6stat[status=="ARI", PIL_6]
  tc.vec <- perIL_6stat[status=="HC2", PIL_6]
  stat.test <- wilcox.test(gt.vec, tc.vec, paired=FALSE, alternative="greater")
  
  ggplot(perIL_6stat, aes(status, PIL_6, fill=status))  + 
    geom_boxplot(width=0.5, color="black", position="dodge") + 
    geom_jitter(aes(status, PIL_6), size=8, height=0, width=0.03) +
    geom_signif(comparisons=list(c("ARI", "HC2")), 
               annotations=(paste("p=", round(stat.test$p.value, digits=4))), 
               map_signif_level=TRUE, size=0.6, textsize=8, tip_length=0, vjust=0.2) + 
    ylab("IL_6+ Naive B cells (%)") +
    ggtitle("IL_6+ Naive B cells") + 
    scale_fill_manual(values=sc, name="Status") +
    theme_bw() +
    theme(text=element_text(size=42), 
          title=element_text(size=40),
          legend.position='none')
  ggsave(paste0(pa, "Perc_IL_6_naiveb_status_stats.png"), width=8, height=12)
  
  #RANKL+ naive B cells
  perRANKLstat <- table(dt[meta=="Naive", .(donor, RANKL_pos, status)]) %>%
    as.data.table() %>%
    .[RANKL_pos==T, NRANKL:=N, by=donor] %>%
    .[, PRANKL:=100*NRANKL/sum(N), by=donor] %>%
    .[, c("N", "NRANKL"):=NULL] %>%
    .[RANKL_pos==T] %>%
    .[, RANKL_pos:=NULL] %>% 
    .[PRANKL!=0.000000] %>%
    .[, donor:=factor(donor, levels=donors)] %>%
    .[, status:=factor(status, levels=clin.status)]
  
  gt.vec <- perRANKLstat[status=="ARI", PRANKL]
  tc.vec <- perRANKLstat[status=="HC2", PRANKL]
  stat.test <- wilcox.test(gt.vec, tc.vec, paired=FALSE, alternative="greater")
  
  ggplot(perRANKLstat, aes(status, PRANKL, fill=status))  + 
    geom_boxplot(width=0.5, color="black", position="dodge") + 
    geom_jitter(aes(status, PRANKL), size=8, height=0, width=0.03) +
    geom_signif(comparisons=list(c("ARI", "HC2")), 
                annotations=(paste("p=", round(stat.test$p.value, digits=4))), 
                map_signif_level=TRUE, size=0.6, textsize=8, tip_length=0, vjust=0.2) + 
    ylab("RANKL+ Naive B cells (%)") +
    ggtitle("RANKL+ Naive B cells") + 
    scale_fill_manual(values=sc, name="Status") +
    theme_bw() +
    theme(text=element_text(size=42), 
          title=element_text(size=40),
          legend.position='none')
  ggsave(paste0(pa, "Perc_RANKL_naiveb_status_stats.png"), width=8, height=12)
  
}


###### MAIN ######

# upload processed and metaclustered B cell cytometry data
dat <- fread(dat.path) %>%
  as.data.table() %>%
  .[, Lin:=NULL] %>%
  .[, Viability:=NULL] %>%
  .[status=="at-risk", status:="ARI"] %>% 
  .[status=="control", status:="HC2"] %>% 
  .[, `:=`(batch=factor(batch, levels=batches), 
           donor=factor(donor, levels=donors),
           cell_type=factor(cell_type), 
           condition=factor(condition),
           time=factor(time), 
           isotype=factor(isotype, levels=isotypes), 
           status=factor(status, levels=clin.status), 
           meta=factor(meta, levels=subsets))]

#Main Figure3I analysis plots
naiveBCytokineStatus(dt=dat)
