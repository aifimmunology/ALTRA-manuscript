####################################################################################################################
#
# Script: IC_flow_bcell_data_process_final.R
# Project: He, Glass et al pre-clinical RA study 
# Subproject: Peripheral B cell intracellular flow cytometry data analysis
# Author: Marla Glass
# Date: 09-04-2025
#
# Process fcs files for analysis of ARI (ACPA+ at-risk individuals) and HC (ACPA- control) B cell flow cytometry data 
#   Input:
#     FCS files gated as viable B Cells
#     CSV file with metadata for fcs files, includes condition and donor/subject ids
#   Output:
#     CSV file with all data and metadata merged
#
# Flow cytometry data pre-processing:
#   Flow PBMC data was unmixed and compensated, then gated on singlets > viable cells > B cells, 
#   Downloaded gated B cell fcs files for each culture condition and sample donor
#
# Data processing in R:
#   Download and annotate live B cells with condition 
#   asinh transform by batch
#   Scale to 99.9th percentile by batch
# 
#######################################################################################################################


###### LIBRARIES ######

require(flowCore)
require(MetaCyto)
require(tidyverse)
require(data.table)
require(listr)

##### INPUTS #####

# *update these file paths as needed*
fcs.path <- "~/data_analysis/fcs/"
table.path <- "~/data_analysis/tables/"

factors <- c("cell_type", "condition", "donor", "time", "sample_status", "batch", "status")

dump <- c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-B-A", "SSC-B-H", "SSC-B-W", "SSC-H", "SSC-W", 
          "Time", "AF-A")

conds <- c("CpG_CD40_stimulated", "unstimulated")
donors <- c("ARI_1" , "ARI_2" , "ARI_3" , "ARI_4" , "ARI_5" , "ARI_6" , 
            "ARI_7" , "ARI_8" , "ARI_9" , "ARI_10" , "ARI_11" , "ARI_12" , 
            "ARI_13" , "ARI_14" , "ARI_15" , "ARI_16" , "ARI_17" , 
            "HC_1" , "HC_2" , "HC_3" , "HC_4" , "HC_5" , "HC_6" , 
            "HC_7" , "HC_8" , "HC_9" , "HC_10" , "HC_11" , "HC_12" , 
            "HC_13", "D_1", "D_2" , "D_3" , "D_4" , "D_5" , "D_6" , 
            "D_7" , "D_8" , "D_9" , "D_10" , "D_11" , "D_12" , 
            "D_13" , "D_14" , "D_15")
sample.statuses <- c("experimental")
batches <- c("1", "2", "3")
clin.status <- c("at-risk", "control")

intracellular <- c("TNFa", "IL_6", "IL_10", "RANKL")
intra.factors <- c("IL_10_pos", "IL_6_pos", "TNFa_pos", "RANKL_pos")

all.exp.data <- c("FSC-A", "FSC-H",	"FSC-W", "SSC-A", "SSC-B-A", "SSC-B-H", "SSC-B-W", "SSC-H", "SSC-W", 
                  "CD107a",	"CD69",	"Lin",	"CD45RB",	"CD24", "CD20",	"HLA_DR",	"CD27",	"CD10",	
                  "CD21", "CD73", "CD95", "CD11c",
                  "IgG", "IgA", "IgD", "IgM", 
                  "TNFa", "IL_6", "IL_10", "RANKL", 
                  "Viability", "Time")

surface.markers <- c("CD107a",	"CD69",	"Lin", "CD45RB", "CD24", "CD20", "HLA_DR", "CD27", "CD10",	
                     "CD21", "CD73", "CD95", "CD11c",
                     "IgG", "IgA", "IgD", "IgM")

##### FUNCTIONS #####

readFiles <- function(p) {
  # takes in a path with fcs files and returns a list of data.tables of the expression matrix
  # Inputs:
  #   p - directory storing fcs files
  # Outputs:
  #   frames - a list of data.tables
  
  print("Reading files")
  
  files <- list.files(path=p, pattern=".fcs", full.names=FALSE, recursive=FALSE)
  frames <- setNames(vector("list", length(files)), files)
  for (i in seq(files)) {
    fcs <- read.FCS(paste0(p, files[i]), transformtablation = FALSE, emptyValue = FALSE)
    frames[[i]] <- data.table(exprs(fcs)) %>%
      setnames(pData(parameters(fcs))$desc)
  }
  
  return(frames)
}

combineFiles1 <- function(frames, fd, du=paste(dump, "CD14"), fa=factors) {
  # Combines a list of data.tables into a single data.table with factor columns added
  # Inputs:
  #   frames - a list of data.tables
  #   fd - csv file location with fcs data
  #   dump - vector of dump channel names
  #   fa - vector of factors
  # Outputs:
  #   frame - a single data.table
  
  print("Combining files")
  
  d <- fread(fd, stringsAsFactors=T)
  all.cols <- unique(unlist(lapply(frames, colnames))) %>% .[!. %in% du] %>% c(fa)
  frame <- data.table(matrix(nrow=0, ncol=length(all.cols))) %>% setnames(all.cols)
  for (f in d$filename) {
    frame <- frames[[f]] %>%
      cbind(d[filename==f]) %>%
      .[, all.cols, with=F] %>%
      rbind(frame)
  }
  
  frame[, `:=`(donor=factor(donor),
               cell_type=factor(cell_type), 
               condition=factor(condition, levels=conds), 
               time=factor(time), 
               batch=factor(batch),
               status=factor(status, levels=clin.status),
               sample_status=factor(sample_status, levels=sample.statuses))] %>%
    setnames(., "IL-6", "IL_6") %>%
    setnames(., "IL-10", "IL_10") %>%
    setnames(., "HLA-DR", "HLA_DR")
  
  return(frame)
}

combineFiles2a3 <- function(frames, fd, du=paste(dump, "CD307d"), fa=factors) {
  # Combines a list of data.tables into a single data.table with factor columns added
  # Inputs:
  #   frames - a list of data.tables
  #   fd - csv file location with fcs data
  #   dump - vector of dump channel names
  #   fa - vector of factors
  # Outputs:
  #   frame - a single data.table
  
  print("Combining files")
  
  d <- fread(fd, stringsAsFactors=T)
  all.cols <- unique(unlist(lapply(frames, colnames))) %>% .[!. %in% du] %>% c(fa)
  frame <- data.table(matrix(nrow=0, ncol=length(all.cols))) %>% setnames(all.cols)
  for (f in d$filename) {
    frame <- frames[[f]] %>%
      cbind(d[filename==f]) %>%
      .[, all.cols, with=F] %>%
      rbind(frame)
  }
  
  frame[, `:=`(donor=factor(donor),
               cell_type=factor(cell_type), 
               condition=factor(condition, levels=conds), 
               time=factor(time), 
               batch=factor(batch),
               status=factor(status, levels=clin.status),
               sample_status=factor(sample_status, levels=sample.statuses))] %>%
    setnames(., "IL6", "IL_6") %>%
    setnames(., "IL10", "IL_10") %>%
    setnames(., "HLADR", "HLA_DR")
  
  return(frame)
}

asinTransform <- function(dt, fa=factors) {
  # asinh data transforms
  # Inputs:
  #   dt - data.table
  #   fa - character vector of channels not to transform
  # Outputs:
  #   dat - data.table
  
  print("Asinh transforming")
  
  to.transform <- colnames(dt)[!colnames(dt) %in% fa]
  dt[, (to.transform) := asinh(dt[, to.transform, with=F]/220),]
  
  return(dt)
}

scaleData <- function(dt, fa=factors, quant.lo=0.001, quant.hi=0.999) {
  # scales data
  # Inputs:
  #   dt - data.table
  #   fa - character vector of channels not to scale
  #   quant.lo - percentile of expression to shift to zero
  #   quant.hi - percentile of expression to scale to
  # Outputs:
  #   dt - data.table
  
  print("Scaling data")
  
  channels <- colnames(dt)[!colnames(dt) %in% fa]
  
  refs.lo <- mapply(quantile, x=dt[, channels, with=F], MARGIN=2, probs=c(quant.lo), na.rm=T)
  refs.hi <- mapply(quantile, x=dt[, channels, with=F], MARGIN=2, probs=c(quant.hi), na.rm=T)
  for (i in seq(channels)) {
    dt[, channels[i]:=(dt[, channels[i], with=F] - refs.lo[i])/ (refs.hi[i] - refs.lo[i])]
  }
  
  return(dt)
}

dataMerge <- function() {
  # Merges flow datasets 
  # Outputs: merged data table
  
  print("merging experiment data tables")
  
  dat.1 <- as.data.table(dat.1) %>%
    .[, CD14:=NULL] %>%
    .[, Time:=NULL] %>%
    .[, `:=`(donor=factor(donor),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             batch=factor(batch),
             status=factor(status, levels=clin.status),
             sample_status=factor(sample_status))]
  dat.2 <- as.data.table(dat.2) %>%
    .[, CD307d:=NULL] %>%
    .[, Time:=NULL] %>%
    .[, `:=`(donor=factor(donor),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             batch=factor(batch),
             status=factor(status, levels=clin.status),
             sample_status=factor(sample_status))]
  dat.3 <- as.data.table(dat.3) %>%
    .[, CD307d:=NULL] %>%
    .[, Time:=NULL] %>%
    .[, `:=`(donor=factor(donor),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             batch=factor(batch),
             status=factor(status, levels=clin.status),
             sample_status=factor(sample_status))]
  
  temp <- merge(dat.1, dat.2, all=T) %>%
    as.data.table() 
  dt <- merge(temp, dat.3, all=T) %>%
    as.data.table() %>%
    .[, `:=`(donor=factor(donor, levels=donors),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             batch=factor(batch, levels=batches),
             status=factor(status, levels=clin.status),
             sample_status=factor(sample_status, levels=sample.statuses))]
  
  return(dt)
}

###### MAIN ######

dat.1 <- readFiles(p=paste0(fcs.path, "batch1/")) %>%
  combineFiles1(fd=paste0(table.path, "ARI_HC_Bcell_process_data_batch1.csv")) %>%
  asinTransform() %>%
  scaleData()
dat.2 <- readFiles(p=paste0(fcs.path, "batch2/")) %>%
  combineFiles2a3(fd=paste0(table.path, "ARI_HC_Bcell_process_data_batch2.csv")) %>%
  asinTransform() %>%
  scaleData() 
dat.3 <- readFiles(p=paste0(fcs.path, "batch3/")) %>%
  combineFiles2a3(fd=paste0(table.path, "ARI_HC_Bcell_process_data_batch3.csv")) %>%
  asinTransform() %>%
  scaleData() 
dat <- dataMerge()
dat[, .N, by=donor]

#save processed single-cell data file
write_csv(dat, paste0(table.path, "bcell_processed_flow_data.csv"))

#optional check of data file
#check.dat <- fread(paste0(table.path, "bcell_processed_flow_data.csv"))
