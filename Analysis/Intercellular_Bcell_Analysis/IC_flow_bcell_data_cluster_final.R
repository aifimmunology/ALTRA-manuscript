##############################################################################
#
# Script: IC_flow_bcell_data_cluster_final.R
# Project: He, Glass et al pre-clinical RA study 
# Subproject: Peripheral B cell intracellular flow cytometry data analysis
# Author: Marla Glass
# Date: 09-04-2025
#
# This program takes in a CSV file of processed viable B cell flow cytometry data from IC_flow_bcell_data_process.R
# Processed flow cytometry data has been arsinh transformed and scaled 
#
# Analyses:
# 1) Label B cells isotypes and set cytokine and RANKL positivity cutoffs, 
# 2) generate clusters and assign metaclusters
#
# Output: CSV files of cytometry data with B cell metaclusters and BCR isotypes assigned  
#
##############################################################################


###### LIBRARIES ######

require(tidyverse)
require(ggplot2)
require(data.table)
require(FlowSOM)
require(viridis)
require(pheatmap)
require(uwot)
require(calecopal)

###### USER INPUTS ######

# *update these file paths as needed*
images.path <- "~/data_analysis/images/"
processimages.path <- "~/data_analysis/images/process/" 
dat.path <- "~/data_analysis/tables/bcell_processed_flow_data.csv"
path <- "~/data_analysis/tables/"

dump <- c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-B-A", "SSC-B-H", "SSC-H", "SSC-W")

factors <- c("cell_type", "condition", "donor", "time", "sample_status", "isotype", "cluster", "meta", "batch", "status", "isotype")

conds <- c("CpG_CD40_stimulated", "unstimulated")

in.donors <- c("ARI_1" , "ARI_2" , "ARI_3" , "ARI_4" , "ARI_5" , "ARI_6" , 
                 "ARI_7" , "ARI_8" , "ARI_9" , "ARI_10" , "ARI_11" , "ARI_12" , 
                 "ARI_13" , "ARI_14" , "ARI_15" , "ARI_16" , "ARI_17" , 
                 "HC_1" , "HC_2" , "HC_3" , "HC_4" , "HC_5" , "HC_6" , 
                 "HC_7" , "HC_8" , "HC_9" , "HC_10" , "HC_11" , "HC_12" , 
                 "HC_13", "D_1", "D_2" , "D_3" , "D_4" , "D_5" , "D_6" , 
                 "D_7" , "D_8" , "D_9" , "D_10" , "D_11" , "D_12" , 
                 "D_13" , "D_14" , "D_15")
donors <- c("ARI_1" , "ARI_2" , "ARI_3" , "ARI_4" , "ARI_5" , "ARI_6" , 
            "ARI_7" , "ARI_8" , "ARI_9" , "ARI_10" , "ARI_11" , "ARI_12" , 
            "ARI_13" , "ARI_14" , "ARI_15" , "ARI_16" , "ARI_17" , 
            "HC_1" , "HC_2" , "HC_3" , "HC_4" , "HC_5" , "HC_6" , 
            "HC_7" , "HC_8" , "HC_9" , "HC_10" , "HC_11" , "HC_12" , 
            "HC_13")

batches <- c("1", "2", "3")

clin.status <- c("at-risk", "control")
clin.colors <- c("#B10906", "#5284A3") %>% setNames(clin.status) 

intracellular <- c("TNFa", "IL_6", "IL_10", "RANKL")
intra.factors <- c("IL_10_pos", "IL_6_pos", "TNFa_pos", "RANKL_pos")

all.cols <- c("CD107a",	"CD69",	"Lin",	"CD45RB",	
              "CD24", "CD20",	"HLA_DR",	"CD27",	"CD10",	
              "CD11c", "CD73", "CD95", "CD11c", "CD21",
              "IgG", "IgA", "IgD", "IgM", 
              "TNFa", "IL_6", "IL_10", "RANKL", 
              "Viability")

surface.markers <- c("CD107a", "CD69", "Lin", "CD45RB",	
                     "CD24", "CD20", "HLA_DR", "CD27", "CD10",	
                     "CD11c", "CD73", "CD95", "CD21",
                     "IgG", "IgA", "IgD", "IgM")

pal <- cal_palette(name="lupinus", n=24, type="continuous")
pal2 <- cal_palette(name="canary", n=21, type="continuous")
donor.colors <- c(pal, pal2) %>% setNames(donors)

subsets <- c("CD27_neg_Effector", "CD27_pos_Effector", "Early_Memory", "Core_Memory", "Late_Memory", "CD95_Memory", 
             "Transitional", "Naive", "Plasma")
subset.colors <- c("#F9B5AC", "#861D31", "#C4E2E1", "#399390", "#849324", "#263FA6",
                   "#805D93", "#CF7C63", "#323949") %>% setNames(subsets)

isotypes <- c("IgD", "IgMD", "IgM", "IgG", "IgA", "ND", "surface_Ig-")
isotype.colors <- c("#664A5B", "#C05746", "#4D7184", "#8DB979", "#012A36", "#D1D2D4", "#ED7D3B") %>% setNames(isotypes)

cluster.colors <- cal_palette("figmtn", 100, type="continuous") %>% as.vector() %>% setNames(1:100)

###### FUNCTIONS ######

processBdata <- function(dt=dat, 
                         markers=surface.markers, 
                         pa=processimages.path) {
  # Further isolate live B cells in flow dataset by removing contaminating PBMC data
  # Inputs:
  #   dt - data.table
  #   subset.markers - character vector of numeric column names; B cell markers only
  #   pa - path to images folder
  # Outputs:
  #   dt - cleaned data.table of live b cell data only 
  
  ### B cell isolation ###
  
  range(dt$CD20)
  range(dt$Lin)
  
  n.sample <- 8000
  set.seed(888)
  subsampled <- dt[, .SD[sample(.N, n.sample)], by=.(batch)]
  
  ggplot(subsampled, aes(CD27, CD20, fill=IgD)) + 
    geom_point(color="black", pch=21) +
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.56)
  
  ggplot(subsampled, aes(CD27, Lin, fill=CD20)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.7)
  
  dt[, non_B:=F]  
  dt[Lin>1, non_B:=T]
  dt[CD20<0.56 & CD27<0.7, non_B:=T]
  dt <- dt[order(non_B)]
  
  n.sample <- 2000
  set.seed(888)
  subsampled <- dt[, .SD[sample(.N, n.sample)], by=.(batch)]
  
  ggplot(subsampled, aes(CD27, CD20, fill=non_B)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + theme()
  ggsave(paste0(pa, "non_B.png"))
  
  #remove non-B cells in the data 
  dt <- dt[non_B==F]
  dt[, non_B:=NULL]
  
  # remove the extra (non-marker) channels in data
  dt[, (dump):=NULL]
  
  dt <- dt[, `:=`(donor=factor(donor, levels=in.donors),
                  cell_type=factor(cell_type), 
                  condition=factor(condition, levels=conds), 
                  time=factor(time), 
                  batch=factor(batch, levels=batches),
                  sample_status=factor(sample_status))]
  
  return(dt)
}
  
isotypeBdata <- function(dt, 
                        markers=surface.markers, 
                        pa=processimages.path) {
    # Labels B cell isotype by experiment
    # Inputs:
    #   dt - data.table
    #   subset.markers - character vector of numeric column names; b cell markers only
    #   pa - path to images folder
    # Outputs:
    #   dt - data.table of live b cell data with added isotype column
    
  ### Isotype labeling - heavy chain ###
  
  # Batch1 isotypes
  min(table(dt$donor))
  n.sample <- 300
  set.seed(888)
  subsampled <- dt[batch=="1", .SD[sample(.N, n.sample)], by=.(donor)]
  
  ggplot(subsampled, aes(CD20, Lin, fill=CD27)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.65)
  
  ggplot(subsampled, aes(IgD, IgM, fill=IgG)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    ylim(-0.1, 1.1) +
    geom_vline(xintercept=0.28) + geom_hline(yintercept=0.33)
  
  ggplot(subsampled, aes(IgD, IgG, fill=IgA)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.28) + geom_hline(yintercept=0.7)
  
  ggplot(subsampled, aes(IgG, IgA, fill=IgM)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.76)
  
  ggplot(subsampled, aes(IgA, IgM, fill=IgD)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.76) + geom_hline(yintercept=0.33)
  
  ggplot(subsampled, aes(IgG, IgM, fill=IgD)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.33)
  
  dt[batch=="1", isotype:="ND"] 
  dt[batch=="1" & IgD>0.28, isotype:="IgD"]
  dt[batch=="1" & IgM>0.33, isotype:="IgM"]
  dt[batch=="1" & (IgM>0.33 & IgD>0.28), isotype:="IgMD"]
  dt[batch=="1" & (IgG>0.7 & IgD<0.28), isotype:="IgG"]
  dt[batch=="1" & (IgA>0.76 & IgD<0.28), isotype:="IgA"]
  dt[batch=="1" & ((IgG<0.7 & IgA<0.76) & (IgD<0.28 & IgM<0.33)), isotype:="surface_Ig-"]
  dt[batch=="1" & (IgG>0.7 & IgA>0.76), isotype:="ND"]
  dt[batch=="1" & (IgG>0.7 & IgD>0.28), isotype:="ND"]
  dt[batch=="1" & (IgG>0.7 & IgM>0.33), isotype:="ND"]
  dt[batch=="1" & (IgA>0.76 & IgD>0.28), isotype:="ND"]
  dt[batch=="1" & (IgA>0.76 & IgM>0.33), isotype:="ND"]
  
  dt <- dt[, `:=` (donor=factor(donor),
                   cell_type=factor(cell_type), 
                   condition=factor(condition, levels=conds), 
                   time=factor(time), 
                   batch=factor(batch),
                   sample_status=factor(sample_status),
                   isotype=factor(isotype, levels=isotypes))]
    #dt[, .N, by=.(isotype)]
  
  
  # Batch2 isotypes
  min(table(dt$donor))
  n.sample <- 300
  set.seed(888)
  subsampled <- dt[batch=="2", .SD[sample(.N, n.sample)], by=.(donor)]
  
  ggplot(subsampled, aes(CD20, Lin, fill=CD27)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.65)
  
  ggplot(subsampled, aes(IgD, IgM, fill=IgA)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.3) + geom_hline(yintercept=0.34)
  
  ggplot(subsampled, aes(IgD, IgG, fill=IgA)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.3) + geom_hline(yintercept=0.77)
  
  ggplot(subsampled, aes(IgG, IgA, fill=IgM)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.77) + geom_hline(yintercept=0.74)
  
  ggplot(subsampled, aes(IgA, IgM, fill=IgD)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.74) + geom_hline(yintercept=0.34)
  
  ggplot(subsampled, aes(IgG, IgM, fill=IgD)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.77) + geom_hline(yintercept=0.34)
  
  dt[batch=="2", isotype:="ND"] 
  dt[batch=="2" & IgD>0.3, isotype:="IgD"]
  dt[batch=="2" & IgM>0.34, isotype:="IgM"]
  dt[batch=="2" & (IgM>0.34 & IgD>0.3), isotype:="IgMD"]
  dt[batch=="2" & (IgG>0.77 & IgD<0.3), isotype:="IgG"]
  dt[batch=="2" & (IgA>0.74 & IgD<0.3), isotype:="IgA"]
  dt[batch=="2" & ((IgG<0.77 & IgA<0.74) & (IgD<0.3 & IgM<0.34)), isotype:="surface_Ig-"]
  dt[batch=="2" & (IgG>0.77 & IgA>0.74), isotype:="ND"]
  dt[batch=="2" & (IgG>0.77 & IgD>0.3), isotype:="ND"]
  dt[batch=="2" & (IgG>0.77 & IgM>0.34), isotype:="ND"]
  dt[batch=="2" & (IgA>0.74 & IgD>0.3), isotype:="ND"]
  dt[batch=="2" & (IgA>0.74 & IgM>0.34), isotype:="ND"]
  
  dt <- dt[, `:=` (donor=factor(donor),
                   cell_type=factor(cell_type), 
                   condition=factor(condition, levels=conds), 
                   time=factor(time), 
                   batch=factor(batch),
                   sample_status=factor(sample_status),
                   isotype=factor(isotype, levels=isotypes))]
    #dt[, .N, by=.(isotype)]
  
  
  # Batch3 isotypes
  n.sample <- 300
  set.seed(888)
  subsampled <- dt[batch=="3", .SD[sample(.N, n.sample)], by=.(donor)]
  
  ggplot(subsampled, aes(CD20, Lin, fill=CD27)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.65)
  
  ggplot(subsampled, aes(IgD, IgM, fill=IgA)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    ylim(0.25, 1.1) +
    xlim(0.4, 1.1) +
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.6)
  
  ggplot(subsampled, aes(IgD, IgG, fill=IgA)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.88)
  
  ggplot(subsampled, aes(IgG, IgA, fill=IgM)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.88) + geom_hline(yintercept=0.78)
  
  ggplot(subsampled, aes(IgA, IgM, fill=IgD)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.78) + geom_hline(yintercept=0.6)
  
  ggplot(subsampled, aes(IgG, IgD, fill=IgD)) + 
    geom_point(color="black", pch=21) + 
    theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.88) + geom_hline(yintercept=0.7)
  
  dt[batch=="3", isotype:="ND"] 
  dt[batch=="3" & IgD>0.7, isotype:="IgD"]
  dt[batch=="3" & IgM>0.6, isotype:="IgM"]
  dt[batch=="3" & (IgM>0.6 & IgD>0.7), isotype:="IgMD"]
  dt[batch=="3" & (IgG>0.88 & IgD<0.7), isotype:="IgG"]
  dt[batch=="3" & (IgA>0.78 & IgD<0.7), isotype:="IgA"]
  dt[batch=="3" & ((IgG<0.88 & IgA<0.78) & (IgD<0.7 & IgM<0.6)), isotype:="surface_Ig-"]
  dt[batch=="3" & (IgG>0.88 & IgA>0.78), isotype:="ND"]
  dt[batch=="3" & (IgG>0.88 & IgD>0.7), isotype:="ND"]
  dt[batch=="3" & (IgG>0.88 & IgM>0.6), isotype:="ND"]
  dt[batch=="3" & (IgA>0.78 & IgD>0.7), isotype:="ND"]
  dt[batch=="3" & (IgA>0.78 & IgM>0.6), isotype:="ND"]
  
  dt <- dt[, `:=` (donor=factor(donor),
                   cell_type=factor(cell_type), 
                   condition=factor(condition, levels=conds), 
                   time=factor(time), 
                   batch=factor(batch),
                   sample_status=factor(sample_status),
                   isotype=factor(isotype, levels=isotypes))]
    #dt[, .N, by=.(isotype)]
  
  
  ## final data
  
  n.sample <- 500
  set.seed(888)
  subsampled <- dt[, .SD[sample(.N, n.sample)], by=.(isotype)]
  
  ggplot(subsampled, aes(CD27, CD45RB, fill=isotype)) + 
    geom_point(color="black", pch=21, size=2) +  
    xlim(0.4, 1.1) + 
    ylim(0, 1.1) + 
    scale_fill_manual(values=isotype.colors) +
    theme_minimal()
  ggsave(paste0(pa, "isotype_b.png"))
  
  dt <- dt[!is.na(isotype)] %>% 
    .[, `:=`(donor=factor(donor, levels=in.donors),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             batch=factor(batch, levels=batches),
             sample_status=factor(sample_status), 
             isotype=factor(isotype, levels=isotypes))]
  
  return(dt)
  
  #save data with B cell isotype assignments as csv file
  #fwrite(dat.b, file=paste0(path, "bcell_processed_flow_data_isotypes.csv"))
}

somCluster <- function(dt, channels, ...) {
  # Clusters with a SOM
  # Inputs:
  #   dt - data.table with all of the data to be clustered
  #   channels - vector of channel names to cluster
  # Outputs:
  #   cluster assignment column
  set.seed(888)
  som.out <- SOM(as.matrix(dt[, channels, with=F]), ...)
  return(as.factor(som.out$mapping[,1]))
}

clusterBData3 <- function(dt=dat[batch=="3"], 
                         subset.markers=setdiff(colnames(dat), c(factors, intracellular, "Lin", "Viability")), 
                         all.markers=setdiff(colnames(dat), c(factors, intracellular, "Viability")), 
                         pa=images.path) {
  # Calls somCluster and then metaclusters flow data into 8 B cell populations for Batch3
  # "CD27_neg_Effector", "CD27_pos_Effector", "Early_Memory", "Core_Memory", "CD95_Memory", "Transitional", "Naive", "Plasma"
  # Inputs:
  #   dt - data.table
  #   subset.markers - character vector of numeric column names
  #   all.markers - character vector of numeric column names
  #   pa - path to images folder
  # Outputs:
  #   dt - data.table with added meta.cluster column
  
  ###
  # Clustering and looking at all clusters
  dt[, cluster:=somCluster(dt, channels=subset.markers, xdim=10, ydim=10)]
  
  medians <- dt[, lapply(.SD, median), .SDcols=all.markers, by=cluster]
  medians.mat <- as.matrix(medians, rownames="cluster") %>% t()
  medians.mat[medians.mat>1] <- 1
  
  #pheatmap(mat=medians.mat, color=magma(20), legend=T, border_color=NA,
  #         filename=paste0(pa, "ARI_HC_3_cluster_heatmap.png"), width=18, height=10)
  #dev.off()
  
  ###
  
  ###
  
  #  Metaclustering
  for (clust in medians$cluster) medians[cluster==clust, count:=nrow(dt[cluster==clust])]
  medians[, adjusted.counts:=log2(count)]
  
  # check and mark remaining CD20- population/cluster 
  ggplot(medians, aes(CD20, Lin, fill=CD27, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + geom_vline(xintercept=0.2) + geom_hline(yintercept=0.88)
  
  ggplot(medians, aes(Lin, HLA_DR, fill=CD20, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + geom_vline(xintercept=0.88) + geom_hline(yintercept=0.5)
  
  medians[Lin>0.88, meta:="non-B"]
  
  # CD27+/- CD11c Effector B cells
  
  ggplot(medians[is.na(meta)], aes(CD11c, CD20, fill=CD21, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.8)
  
  ggplot(medians[is.na(meta)], aes(CD11c, CD21, fill=CD20, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.5)
  
  ggplot(medians[is.na(meta)], aes(CD11c, CD27, fill=CD20, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.68)

  medians[is.na(meta) & (CD11c>0.71 & CD20>0.8 & CD21<0.5 & CD27>0.68), meta:="CD27_pos_Effector"]
  medians[is.na(meta) & (CD11c>0.71 & CD20>0.8 & CD21<0.5 & CD27<0.68), meta:="CD27_neg_Effector"]
  
  # Memory (experienced) and non-memory (inexperienced) B cells
  
  ggplot(medians[is.na(meta)], aes(CD27, IgD, fill=IgM, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.66) + geom_hline(yintercept=0.7)
  
  ggplot(medians[is.na(meta)], aes(CD27, IgA, fill=CD45RB, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.66) + geom_hline(yintercept=0.76)
  
  ggplot(medians[is.na(meta)], aes(CD27, IgG, fill=IgD, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.66) + geom_hline(yintercept=0.7)
  
  ggplot(medians[is.na(meta)], aes(CD27, IgM, fill=IgD, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.66) + geom_hline(yintercept=0.6)
  
  ggplot(medians[is.na(meta)], aes(CD45RB, IgM, fill=IgA, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.74) + geom_hline(yintercept=0.6)
  
  ggplot(medians[is.na(meta)], aes(CD45RB, IgA, fill=CD27, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.74) + geom_hline(yintercept=0.76)
  
  ggplot(medians[is.na(meta)], aes(CD27, CD45RB, fill=CD95, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    xlim(0.5, 1.1) + ylim(0.5, 1.1) +
    geom_vline(xintercept=0.66) + geom_hline(yintercept=0.74)
  
  medians[is.na(meta) & (CD27>0.66 | IgA>0.76 | CD45RB>0.74), meta:="Memory"]
  medians[is.na(meta), meta:="Inexperienced"]
  
  # Transitional (inexperienced) B cell subset
  
  ggplot(medians[meta=="Inexperienced"], aes(CD10, CD24, fill=IgM, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") +
    xlim(0.7,0.85) + ylim(0.5,0.9) +
    scale_size_continuous(range = c(2, 18)) + theme_bw() +  scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.76) + geom_hline(yintercept=0.75)
  
  medians[meta=="Inexperienced" & CD10>0.76 & CD24>0.75, meta:="Transitional"]
  
  # Naive (inexperienced) B cell subsets
  
  ggplot(medians[meta=="Inexperienced"], aes(CD45RB, CD27, fill=IgM, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + 
    scale_fill_viridis(option="B") + 
    theme(legend.position="none") 
  
  ggplot(medians[meta=="Inexperienced"], aes(CD45RB, CD73, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_hline(yintercept=0.7)
  
  medians[meta=="Inexperienced" & CD45RB<0.75, meta:="Naive"]
  
  #Plasma cell subset
  
  ggplot(medians[meta=="Memory"], aes(CD20, CD27, fill=CD45RB, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.62) + geom_hline(yintercept=0.7)
  
  ggplot(medians[meta=="Memory"], aes(CD20, CD27, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.62) + geom_hline(yintercept=0.7)
  
  ggplot(medians[meta=="Memory"], aes(CD20, IgG, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.62) + geom_hline(yintercept=0.74)
  
  medians[meta=="Memory" & CD27>0.7 & CD20<0.62, meta:="Plasma"]
  
  #CD95 Memory B cell subset
  
  ggplot(medians[meta=="Memory"], aes(CD95, CD27, fill=CD45RB, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.80) + geom_hline(yintercept=0.6)
  
  ggplot(medians[meta=="Memory"], aes(CD95, IgG, fill=CD11c, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.8) + geom_hline(yintercept=0.76)
  
  medians[meta=="Memory" & CD27>0.6 & CD95>0.8, meta:="CD95_Memory"]
  
  # CD45RB+ CD27-, CD45RB- CD27+, CD45RB+ CD27+ Memory (experienced) B cell subsets
  
  ggplot(medians[meta=="Memory"], aes(IgG, CD27, fill=IgG, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.76) + geom_hline(yintercept=0.66)
  
  ggplot(medians[meta=="Memory"], aes(CD27, CD45RB, fill=CD95, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.66) + geom_hline(yintercept=0.71) 
   
   ggplot(medians[meta=="Memory"], aes(CD27, CD45RB, fill=IgG, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.66) + geom_hline(yintercept=0.71)
  
  medians[meta=="Memory" & CD27>0.66 & CD45RB<0.71, meta:="Late_Memory"]
  medians[meta=="Memory" & CD27<0.66 & CD45RB>0.71, meta:="Early_Memory"]
  medians[meta=="Memory" & CD27>0.66 & CD45RB>0.71, meta:="Core_Memory"]
  
  medians[meta=="Memory" & IgD<0.5, meta:="Other_Memory"]
  
  #check unlabeled cell clusters
  ggplot(medians[meta=="Inexperienced"], aes(CD20, CD27, fill=IgM, size=count)) + 
    geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.5) + geom_hline(yintercept=0.5)
  
  # medians[meta=="Inexperienced" & CD27>0.4, meta:="Other_Memory"]
  
  table(medians$meta)
  medians[, meta:=factor(meta, levels=subsets)]
  for (clust in medians$cluster) dt[cluster==clust, meta:=medians[cluster==clust, meta]]
  
  # Delete non-B cells from the data.table with metaclusters (using cluster name)
  dt[meta!="Non-B"]
  
  #View expression heatmap with metacluster annotation for unstimulated B cells
  clust.mat <- as.matrix(medians[, !c("meta", "count", "adjusted.counts")], rownames="cluster") %>% t()
  clust.mat[clust.mat>1] <- 1
  clust.annot <- as.data.frame(medians[, .(meta, cluster)])
  rownames(clust.annot) <- clust.annot$cluster
  clust.annot$cluster <- NULL
  ac <- list(meta=subset.colors)
  
  #pheatmap(mat=clust.mat, annotation_col=clust.annot, annotation_colors=ac, 
  #         color=magma(50), legend=F, border_color=NA, 
  #         filename=paste0(pa, "ARI_HC_3_cluster_heatmap_final.png"), width=20, height=10)
  #dev.off()
  
  #final clustered data
  dt3 <- dt[!is.na(meta)] %>% 
    .[, cluster:=NULL] %>%
    .[, `:=`(donor=factor(donor),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             batch=factor(batch), 
             sample_status=factor(sample_status),
             isotype=factor(isotype, levels=isotypes), 
             meta=factor(meta, levels=subsets))]
  #optional: save data with B cell cluster/metacluster assignments as csv file
  #fwrite(dt3, file=paste0(path, "bcell_clustered_flow_data_3.csv"))
}

clusterBData2 <- function(dt=dat[batch=="2"], 
                               subset.markers=setdiff(colnames(dat), c(factors, intracellular, "Lin", "Viability")), 
                               all.markers=setdiff(colnames(dat), c(factors, intracellular, "Viability")), 
                               pa=images.path) {
  # Calls somCluster and then metaclusters flow data into 8 B cell populations for Batch2
  # "CD27_neg_Effector", "CD27_pos_Effector", "Early_Memory", "Core_Memory", "CD95_Memory", "Transitional", "Naive", "Plasma"
  # Inputs:
  #   dt - data.table
  #   subset.markers - character vector of numeric column names
  #   all.markers - character vector of numeric column names
  #   pa - path to images folder
  # Outputs:
  #   dt - data.table with added meta.cluster column

  ###
  # Clustering and looking at all clusters
  dt[, cluster:=somCluster(dt, channels=subset.markers, xdim=10, ydim=10)]
  
  medians <- dt[, lapply(.SD, median), .SDcols=all.markers, by=cluster]
  medians.mat <- as.matrix(medians, rownames="cluster") %>% t()
  medians.mat[medians.mat>1] <- 1
  
  #pheatmap(mat=medians.mat, color=magma(20), legend=T, border_color=NA,
  #         filename=paste0(pa, "ARI_HC_2_cluster_heatmap.png"), width=18, height=10)
  #dev.off()
  
  ###
  
  ###
  
  #  Metaclustering
  for (clust in medians$cluster) medians[cluster==clust, count:=nrow(dt[cluster==clust])]
  medians[, adjusted.counts:=log2(count)]
  
  # check and mark remaining CD20- population/cluster 
  ggplot(medians, aes(CD20, Lin, fill=CD27, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + geom_vline(xintercept=0.2) + geom_hline(yintercept=0.88)
  
  ggplot(medians, aes(Lin, HLA_DR, fill=CD20, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + geom_vline(xintercept=0.88) + geom_hline(yintercept=0.5)
  
  medians[Lin>0.88, meta:="non-B"]
  
  # CD27+/- CD11c Effector B cells
  
  ggplot(medians[is.na(meta)], aes(CD11c, CD20, fill=CD21, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    ylim(0.65, 1) +
    theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.8)
  
  ggplot(medians[is.na(meta)], aes(CD11c, CD21, fill=CD20, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.6)
  
  ggplot(medians[is.na(meta)], aes(CD11c, CD27, fill=CD21, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    ylim(0.5, 1) +
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.69)
  
  ggplot(medians[is.na(meta)], aes(CD21, CD27, fill=CD11c, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.6) + geom_hline(yintercept=0.69)
  
  medians[is.na(meta) & (CD11c>0.7 & CD20>0.8 & CD21<0.6 & CD27>0.69), meta:="CD27_pos_Effector"]
  medians[is.na(meta) & (CD11c>0.7 & CD20>0.8 & CD21<0.6 & CD27<0.69), meta:="CD27_neg_Effector"]
  
  # Memory (experienced) and non-memory (inexperienced) B cells
  
  ggplot(medians[is.na(meta)], aes(CD27, IgD, fill=IgM, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.7)
  
  ggplot(medians[is.na(meta)], aes(CD27, IgA, fill=CD45RB, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.76)
  
  ggplot(medians[is.na(meta)], aes(CD27, IgG, fill=IgD, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.75)
  
  ggplot(medians[is.na(meta)], aes(CD27, IgM, fill=IgA, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.28)
  
  ggplot(medians[is.na(meta)], aes(CD45RB, IgM, fill=IgA, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.75) + geom_hline(yintercept=0.28)
  
  ggplot(medians[is.na(meta)], aes(CD45RB, IgA, fill=CD27, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.8) + geom_hline(yintercept=0.76)
  
  ggplot(medians[is.na(meta)], aes(CD27, CD45RB, fill=CD24, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    xlim(0.5, 1.1) + ylim(0.5, 1.1) +
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.8)
  
  medians[is.na(meta) & (CD27>0.71 | IgA>0.76 | CD45RB>0.8), meta:="Memory"]
  medians[is.na(meta), meta:="Inexperienced"]
  
  # Transitional (inexperienced) B cell subset
  
  ggplot(medians[meta=="Inexperienced"], aes(CD10, CD24, fill=IgM, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") +
   xlim(0.6,0.85) + ylim(0.6,0.95) +
    scale_size_continuous(range = c(2, 18)) + theme_bw() +  scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.8)
  
  medians[meta=="Inexperienced" & CD10>0.71 & CD24>0.8, meta:="Transitional"]
  
  # Naive (inexperienced) B cell subsets
  
  ggplot(medians[meta=="Inexperienced"], aes(CD24, CD27, fill=IgM, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + 
    scale_fill_viridis(option="B") + 
    theme(legend.position="none") + geom_hline(yintercept=0.7)
  
  ggplot(medians[meta=="Inexperienced"], aes(CD24, CD73, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.85) + geom_hline(yintercept=0.76)
  
  medians[meta=="Inexperienced" & CD27<0.7, meta:="Naive"]
  
  #Plasma cell subset
  
  ggplot(medians[meta=="Memory"], aes(CD20, CD27, fill=CD45RB, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.6) + geom_hline(yintercept=0.75)
  
  ggplot(medians[meta=="Memory"], aes(CD20, CD27, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.6) + geom_hline(yintercept=0.75)
  
  ggplot(medians[meta=="Memory"], aes(CD20, IgA, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.6) + geom_hline(yintercept=0.68)
  
  medians[meta=="Memory" & CD27>0.75 & CD20<0.6, meta:="Plasma"]
  
  #CD95 Memory B cell subset
  
  ggplot(medians[meta=="Memory"], aes(CD95, CD27, fill=CD45RB, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.5) + geom_hline(yintercept=0.6)
  
  ggplot(medians, aes(CD95, CD11c, fill=CD21, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.5) + geom_hline(yintercept=0.7)
  
  medians[meta=="Memory" & CD27>0.6 & CD95>0.5, meta:="CD95_Memory"]
  
  # CD45RB+ CD27-, CD45RB- CD27+, CD45RB+ CD27+ Memory (experienced) B cell subsets
  
  ggplot(medians[meta=="Memory"], aes(IgG, CD27, fill=IgM, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.67) + geom_hline(yintercept=0.71)
  
  ggplot(medians[meta=="Memory"], aes(CD27, CD45RB, fill=CD95, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.8) 
  
  ggplot(medians[meta=="Memory"], aes(CD27, CD45RB, fill=IgG, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.71) + geom_hline(yintercept=0.8) 
  
  medians[meta=="Memory" & CD27>0.71 & CD45RB<0.8, meta:="Late_Memory"]
  medians[meta=="Memory" & CD27<0.71 & CD45RB>0.8, meta:="Early_Memory"]
  medians[meta=="Memory" & CD27>0.71 & CD45RB>0.8, meta:="Core_Memory"]
  
  medians[meta=="Memory" & IgD<0.5, meta:="Other_Memory"]
  
  #check unlabeled cell clusters
  ggplot(medians[meta=="Memory"], aes(CD20, CD27, fill=IgM, size=count)) + 
    geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.5) + geom_hline(yintercept=0.5)
  
  # medians[meta=="Inexperienced" & CD27>0.7, meta:="Other_Memory"]
  
  table(medians$meta)
  medians[, meta:=factor(meta, levels=subsets)]
  for (clust in medians$cluster) dt[cluster==clust, meta:=medians[cluster==clust, meta]]
  
  # Delete non-B cells from the data.table with metaclusters (using cluster name)
  dt[meta!="Non-B"]
  
  #View expression heatmap with metacluster annotation for unstimulated B cells
  clust.mat <- as.matrix(medians[, !c("meta", "count", "adjusted.counts")], rownames="cluster") %>% t()
  clust.mat[clust.mat>1] <- 1
  clust.annot <- as.data.frame(medians[, .(meta, cluster)])
  rownames(clust.annot) <- clust.annot$cluster
  clust.annot$cluster <- NULL
  ac <- list(meta=subset.colors)
  
  #pheatmap(mat=clust.mat, annotation_col=clust.annot, annotation_colors=ac, 
  #         color=magma(50), legend=F, border_color=NA, 
  #         filename=paste0(pa, "ARI_HC_2_cluster_heatmap_final.png"), width=20, height=10)
  #dev.off()
  
  #final clustered data
  dt2 <- dt[!is.na(meta)] %>% 
    .[, cluster:=NULL] %>%
    .[, `:=`(donor=factor(donor),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             batch=factor(batch), 
             sample_status=factor(sample_status),
             isotype=factor(isotype, levels=isotypes), 
             meta=factor(meta, levels=subsets))]
  #optional: save data with B cell cluster/metacluster assignments as csv file
  #fwrite(dt2, file=paste0(path, "bcell_clustered_flow_data_2.csv"))
}

clusterBData1 <- function(dt=dat[batch=="1"], 
                               subset.markers=setdiff(colnames(dat), c(factors, intracellular, "Lin", "Viability")), 
                               all.markers=setdiff(colnames(dat), c(factors, intracellular, "Viability")), 
                               pa=images.path) {
  # Calls somCluster and then metaclusters flow data into 8 B cell populations for Batch1
  # "CD27_neg_Effector", "CD27_pos_Effector", "Early_Memory", "Core_Memory", "CD95_Memory", "Transitional", "Naive", "Plasma"
  # Inputs:
  #   dt - data.table
  #   subset.markers - character vector of numeric column names
  #   all.markers - character vector of numeric column names
  #   pa - path to images folder
  # Outputs:
  #   dt - data.table with added meta.cluster column
  
  ###
  # Clustering and looking at all clusters
  dt[, cluster:=somCluster(dt, channels=subset.markers, xdim=10, ydim=10)]
  
  medians <- dt[, lapply(.SD, median), .SDcols=all.markers, by=cluster]
  medians.mat <- as.matrix(medians, rownames="cluster") %>% t()
  medians.mat[medians.mat>1] <- 1
  
  #pheatmap(mat=medians.mat, color=magma(20), legend=T, border_color=NA,
  #         filename=paste0(pa, "ARI_HC_1_cluster_heatmap.png"), width=18, height=10)
  # dev.off()
  
  ###
  
  ###
  
  #  Metaclustering
  for (clust in medians$cluster) medians[cluster==clust, count:=nrow(dt[cluster==clust])]
  medians[, adjusted.counts:=log2(count)]
  
  # check and mark remaining CD20- population/cluster 
  ggplot(medians, aes(CD20, Lin, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + geom_vline(xintercept=0.2) + geom_hline(yintercept=0.88)
  
  ggplot(medians, aes(Lin, HLA_DR, fill=CD20, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + geom_vline(xintercept=0.8) + geom_hline(yintercept=0.3)
  
  medians[Lin>0.8, meta:="non-B"]
  
  # CD27+/- CD11c Effector B cells
  
  ggplot(medians[is.na(meta)], aes(CD11c, CD20, fill=CD21, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") +  ylim(0.65, 1) +
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.81)
  
  ggplot(medians[is.na(meta)], aes(CD11c, CD21, fill=CD20, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.6)
  
  ggplot(medians[is.na(meta)], aes(CD11c, CD27, fill=CD20, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.68)
  
  ggplot(medians[is.na(meta)], aes(CD21, CD27, fill=CD11c, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.5, vjust=0, size=6, color="dark gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.6) + geom_hline(yintercept=0.68)

  medians[is.na(meta) & (CD11c>0.68 & CD20>0.81 & CD21<0.6 & CD27>0.68), meta:="CD27_pos_Effector"]
  medians[is.na(meta) & (CD11c>0.68 & CD20>0.81 & CD21<0.6 & CD27<0.68), meta:="CD27_neg_Effector"]
  
  # Memory (experienced) and non-memory (inexperienced) B cells
  
  ggplot(medians[is.na(meta)], aes(CD27, IgD, fill=IgM, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.6)
  
  ggplot(medians[is.na(meta)], aes(CD27, IgA, fill=CD45RB, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.76)
  
  ggplot(medians[is.na(meta)], aes(CD27, IgG, fill=IgD, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.7)
  
  ggplot(medians[is.na(meta)], aes(CD27, IgM, fill=IgD, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.3)
  
  ggplot(medians[is.na(meta)], aes(CD45RB, IgM, fill=IgA, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.77) + geom_hline(yintercept=0.6)
  
  ggplot(medians[is.na(meta)], aes(CD45RB, IgA, fill=CD27, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.77) + geom_hline(yintercept=0.76)
  
  ggplot(medians[is.na(meta)], aes(CD27, CD45RB, fill=IgD, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    xlim(0.5, 1) + ylim(0.5, 1) +
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.77)
  
  medians[is.na(meta) & (CD27>0.68 | IgA>0.76 | CD45RB>0.77), meta:="Memory"]
  medians[is.na(meta), meta:="Inexperienced"]
  
  # Transitional (inexperienced) B cell subset
  
  ggplot(medians[meta=="Inexperienced"], aes(CD10, CD24, fill=IgM, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") +
    xlim(0.6,0.85) + ylim(0.4,0.9) +
    scale_size_continuous(range = c(2, 18)) + theme_bw() +  scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_vline(xintercept=0.75) + geom_hline(yintercept=0.74)
  
  medians[meta=="Inexperienced" & CD10>0.75 & CD24>0.74, meta:="Transitional"]
  
  # Naive (inexperienced) B cell subsets
  
  ggplot(medians[meta=="Inexperienced"], aes(CD45RB, CD27, fill=IgM, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + 
    scale_fill_viridis(option="B") + 
    theme(legend.position="none") + geom_hline(yintercept=0.7)
  
  ggplot(medians[meta=="Inexperienced"], aes(CD45RB, CD73, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + 
    theme(legend.position="none") + 
    geom_hline(yintercept=0.75)
  
  medians[meta=="Inexperienced" & CD27<0.7, meta:="Naive"]
  
  #Plasma cell subset
  
  ggplot(medians[meta=="Memory"], aes(CD20, CD27, fill=CD45RB, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.7)
  
  ggplot(medians[meta=="Memory"], aes(CD20, CD27, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.7)
  
  ggplot(medians[meta=="Memory"], aes(CD20, IgG, fill=HLA_DR, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.7)
  
  medians[meta=="Memory" & CD27>0.7 & CD20<0.7, meta:="Plasma"]
  
  #CD95 Memory B cell subset
  
  ggplot(medians[meta=="Memory"], aes(CD95, CD27, fill=CD11c, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.54) + geom_hline(yintercept=0.68)
  
  ggplot(medians[meta=="Memory"], aes(CD95, IgG, fill=CD11c, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.54) + geom_hline(yintercept=0.76)
  
  medians[meta=="Memory" & CD27>0.68 & CD95>0.54, meta:="CD95_Memory"]
  
  # CD45RB+ CD27-, CD45RB- CD27+, CD45RB+ CD27+ Memory (experienced) B cell subsets
  
  ggplot(medians[meta=="Memory"], aes(IgG, CD27, fill=IgG, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.7) + geom_hline(yintercept=0.68)
  
  ggplot(medians[meta=="Memory"], aes(CD27, CD45RB, fill=IgM, size=count)) + geom_point(color="black", pch=21) +
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="gray") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.77) 
  
  ggplot(medians[meta=="Memory"], aes(CD27, CD45RB, fill=IgG, size=count)) + geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.68) + geom_hline(yintercept=0.77)
  
  medians[meta=="Memory" & CD27>0.68 & CD45RB<0.77, meta:="Late_Memory"]
  medians[meta=="Memory" & CD27<0.68 & CD45RB>0.77, meta:="Early_Memory"]
  medians[meta=="Memory" & CD27>0.68 & CD45RB>0.77, meta:="Core_Memory"]
  
  medians[meta=="Memory" & IgD<0.5, meta:="Other_Memory"]
  
  #check unlabeled cell clusters
  ggplot(medians[meta=="Inexperienced"], aes(CD20, CD27, fill=IgM, size=count)) + 
    geom_point(color="black", pch=21) + 
    geom_text(aes(label=cluster), hjust=0.45, vjust=0, size=6, color="black") + 
    scale_size_continuous(range = c(2, 18)) + theme_bw() + scale_fill_viridis(option="B") + theme(legend.position="none") + 
    geom_vline(xintercept=0.5) + geom_hline(yintercept=0.5)
  
  # medians[meta=="Inexperienced" & CD27>0.4, meta:="Other_Memory"]
  
  table(medians$meta)
  medians[, meta:=factor(meta, levels=subsets)]
  for (clust in medians$cluster) dt[cluster==clust, meta:=medians[cluster==clust, meta]]
  
  # Delete non-B cells from the data.table with metaclusters (using cluster name)
  dt[meta!="Non-B"]
  
  #View expression heatmap with metacluster annotation for unstimulated B cells
  clust.mat <- as.matrix(medians[, !c("meta", "count", "adjusted.counts")], rownames="cluster") %>% t()
  clust.mat[clust.mat>1] <- 1
  clust.annot <- as.data.frame(medians[, .(meta, cluster)])
  rownames(clust.annot) <- clust.annot$cluster
  clust.annot$cluster <- NULL
  ac <- list(meta=subset.colors)
  
  #pheatmap(mat=clust.mat, annotation_col=clust.annot, annotation_colors=ac, 
  #         color=magma(50), legend=F, border_color=NA, 
  #          filename=paste0(pa, "ARI_HC_1_cluster_heatmap_final.png"), width=20, height=10)
  # dev.off()
  
  #final clustered data
  dt1 <- dt[!is.na(meta)] %>% 
    .[, cluster:=NULL] %>%
    .[, `:=`(donor=factor(donor),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             batch=factor(batch), 
             sample_status=factor(sample_status),
             isotype=factor(isotype, levels=isotypes), 
             meta=factor(meta))]
  #optional: save data with B cell cluster/metacluster assignments as csv file
  #fwrite(dt1, file=paste0(path, "bcell_clustered_flow_data_1.csv"))
}

setCutoffsUnstim <- function(dt=dat, 
                             cyt=intracellular, 
                             quant=0.99,
                             pa=processimages.path) {
  # Identify cutoffs for cytokines and RANKL and labels cells based on percentile
  
  # Inputs:
  #   dt - data.table
  #   cyt - character vector of column names to use
  #   quant - percentile of unstimulated condition to choose as the positivity cutoff
  # Outputs:
  #   dt - data.table
  
  print("Setting cutoffs for cytokine and RANKL positivity")
  
  dt[, `:=`(TNFa_pos=F, IL_6_pos=F, IL_10_pos=F, RANKL_pos=F)]
  
  for (d in in.donors) {
    for (m in subsets) {
      
      cutoffs.per <- dt[meta==m & donor==d & condition=="unstimulated", lapply(.SD, quantile, probs=c(quant), na.rm=T), .SDcols=cyt] %>% 
        as.list() %>% unlist()
      
      dt[donor==d & meta==m & TNFa>cutoffs.per["TNFa"], TNFa_pos:=T]
      dt[donor==d & meta==m & IL_6>cutoffs.per["IL_6"], IL_6_pos:=T]
      dt[donor==d & meta==m & IL_10>cutoffs.per["IL_10"], IL_10_pos:=T]
      dt[donor==d & meta==m & RANKL>cutoffs.per["RANKL"], RANKL_pos:=T]
    }
  }
  
  n.sample <- 3000
  set.seed(888)
  subsampled <- dt[, .SD[sample(.N, n.sample)], by=.(batch)] 
  
  ggplot(subsampled[CD20>0.65 & batch=="2"], aes(CD20, TNFa, fill=TNFa_pos)) + 
    geom_point(color="black", pch=21, size=6) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    xlim(0, 1.0) +
    ylim(0, 1.0) +
    theme_minimal() +
    theme(legend.position="none", 
          text=element_text(size=30),
          axis.title=element_text(size=40))
 #ggsave(paste0(pa, "TNFa_pos_B_cell_biaxial.png"), width=8, height=8)
 
   ggplot(subsampled[condition=="CpG_CD40_stimulated"], aes(CD20, TNFa, fill=TNFa_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  ggsave(paste0(pa, "TNFa_pos_stim_B_cell_biaxial.png"), width=8, height=6)
  ggplot(subsampled[condition=="unstimulated"], aes(CD20, TNFa, fill=TNFa_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
 #ggsave(paste0(pa, "TNFa_pos_unstim_B_cell_biaxial.png"), width=8, height=6)
  
  ggplot(subsampled, aes(CD20, IL_6, fill=IL_6_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  #ggsave(paste0(pa, "IL_6_pos_B_cell_biaxial.png"), width=8, height=6)
  ggplot(subsampled[condition=="CpG_CD40_stimulated"], aes(CD20, IL_6, fill=IL_6_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  #ggsave(paste0(pa, "IL_6_pos_stim_B_cell_biaxial.png"), width=8, height=6)
  ggplot(subsampled[condition=="unstimulated"], aes(CD20, IL_6, fill=IL_6_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  #ggsave(paste0(pa, "IL_6_pos_unstim_B_cell_biaxial.png"), width=8, height=6)
  
  ggplot(subsampled, aes(CD20, IL_10, fill=IL_10_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  #ggsave(paste0(pa, "IL_10_pos_B_cell_biaxial.png"), width=8, height=6)
  ggplot(subsampled[condition=="CpG_CD40_stimulated"], aes(CD20, IL_10, fill=IL_10_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  #ggsave(paste0(pa, "IL_10_pos_stim_B_cell_biaxial.png"), width=8, height=6)
  ggplot(subsampled[condition=="unstimulated"], aes(CD20, IL_10, fill=IL_10_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  #ggsave(paste0(pa, "IL_10_pos_unstim_B_cell_biaxial.png"), width=8, height=6)
  
  ggplot(subsampled, aes(CD20, RANKL, fill=RANKL_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  #ggsave(paste0(pa, "RANKL_pos_B_cell_biaxial.png"), width=8, height=6)
  ggplot(subsampled[condition=="CpG_CD40_stimulated"], aes(CD20, RANKL, fill=RANKL_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  #ggsave(paste0(pa, "RANKL_pos_stim_B_cell_biaxial.png"), width=8, height=6)
  ggplot(subsampled[condition=="unstimulated"], aes(CD20, RANKL, fill=RANKL_pos)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_manual(values=c('#EFBFBD', '#2F4858')) +
    theme(legend.position="right") + 
    theme_minimal() 
  #ggsave(paste0(pa, "RANKL_pos_unstim_B_cell_biaxial.png"), width=8, height=6)
  
  return(dt)
}


###### MAIN ######

#upload processed flow viable B cell data
dat <- fread(dat.path) %>%
  .[, `:=`(donor=factor(donor, levels=in.donors),
           cell_type=factor(cell_type), 
           condition=factor(condition, levels=conds), 
           time=factor(time), 
           batch=factor(batch, levels=batches), 
           sample_status=factor(sample_status))]
# B cell isolation and isotype labeling
dat.b <- processBdata() %>% 
  isotypeBdata()
dat.b[, .N, by=isotype]
#save data with purely B cells and isotype assignments as csv file
fwrite(dat.b, file=paste0(path, "bcell_processed_flow_data_isotypes.csv"))
#proceed with cleaned and isotype-labeled B cell data 
dat <- dat.b %>%
  .[, `:=`(donor=factor(donor, levels=in.donors),
           cell_type=factor(cell_type), 
           condition=factor(condition, levels=conds), 
           time=factor(time), 
           batch=factor(batch, levels=batches), 
           sample_status=factor(sample_status),
           isotype=factor(isotype, levels=isotypes))] 
#cluster flow B cell data
clusterBData3()
clusterBData2()
clusterBData1()
temp <- merge(dt3, dt2, all=T) %>%
  as.data.table() 
dat <- merge(temp, dt1, all=T) %>%
  as.data.table() %>%
  .[, `:=`(donor=factor(donor, levels=in.donors),
           cell_type=factor(cell_type), 
           condition=factor(condition, levels=conds), 
           time=factor(time), 
           batch=factor(batch, levels=batches), 
           sample_status=factor(sample_status),
           isotype=factor(isotype, levels=isotypes), 
           meta=factor(meta, levels=subsets))]
dat[, .N, by=.(batch, meta)]
sum((dt3[,.N]), (dt2[,.N]), (dt1[,.N]))
dat[,.N]
dat <- setCutoffsUnstim() %>% 
  .[status!="na"] %>% 
  .[, `:=`(donor=factor(donor, levels=donors),
           cell_type=factor(cell_type), 
           condition=factor(condition, levels=conds), 
           time=factor(time), 
           batch=factor(batch, levels=batches), 
           status=factor(status, levels=clin.status),
           sample_status=factor(sample_status),
           isotype=factor(isotype, levels=isotypes), 
           meta=factor(meta, levels=subsets))] %>% 
  .[!is.na(donor)]

#save single-cell data file with B cell metacluster assignments
fwrite(dat, file=paste0(path, "bcell_clustered_flow_data.csv"))
