####################################################################################################################
#
# Script: IC_flow_bcell_data_process.R
# Project: He, Glass et al pre-clinical RA study 
# Subproject: Peripheral B cell intracellular flow cytometry data analysis
# Author: Marla Glass
# Date: 07-31-25
#
# Process fcs files for analysis of ARI (ACPA+) and HC2/CON (ACPA- control) B cell flow cytometry data 
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
#   asinh transform by experiment
#   Scale to 99.9th percentile by experiment
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

factors <- c("cell_type", "condition", "donor", "time", "sample_status", "expt")

dump <- c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-B-A", "SSC-B-H", "SSC-B-W", "SSC-H", "SSC-W", 
          "Time", "AF-A")

conds <- c("CpG_CD40_stimulated", "unstimulated")
donors <- c("ARI_1" , "ARI_2" , "ARI_3" , "ARI_4" , "ARI_5" , "ARI_6" , "ARI_7" , "ARI_8" , "ARI_9" , "ARI_10" , "ARI_11" , "ARI_12" , 
            "ARI_13" , "ARI_14" , "ARI_15" , "ARI_16" , "ARI_17" , "ARI_18" , "ARI_19" , "ARI_20" , "ARI_21" , "ARI_22" , "ARI_23" , "ARI_24" , 
            "CON_1" , "CON_2" , "CON_3" , "CON_4" , "CON_5" , "CON_6" , "CON_7" , "CON_8" , "CON_9" , "CON_10" , "CON_11" , "CON_12" , 
            "CON_13" , "CON_14" , "CON_15" , "CON_16" , "CON_17" , "CON_18" , "CON_19" , "CON_20" , "CON_21")
sample.statuses <- c("experimental")
expts <- c("386", "490", "537")

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

combineFiles386 <- function(frames, fd, du=paste(dump, "CD14"), fa=factors) {
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
               expt=factor(expt),
               sample_status=factor(sample_status, levels=sample.statuses))] %>%
    setnames(., "IL-6", "IL_6") %>%
    setnames(., "IL-10", "IL_10") %>%
    setnames(., "HLA-DR", "HLA_DR")
  
  return(frame)
}

combineFiles490a537 <- function(frames, fd, du=paste(dump, "CD307d"), fa=factors) {
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
               expt=factor(expt),
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

exptDataMerge <- function() {
  # Merges flow datasets 
  # Outputs: merged data table
  
  print("merging experiment data tables")
  
  dat.386 <- as.data.table(dat.386) %>%
    .[, CD14:=NULL] %>%
    .[, Time:=NULL] %>%
    .[, `:=`(donor=factor(donor),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             expt=factor(expt),
             sample_status=factor(sample_status))]
  dat.490 <- as.data.table(dat.490) %>%
    .[, CD307d:=NULL] %>%
    .[, Time:=NULL] %>%
    .[, `:=`(donor=factor(donor),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             expt=factor(expt),
             sample_status=factor(sample_status))]
  dat.537 <- as.data.table(dat.537) %>%
    .[, CD307d:=NULL] %>%
    .[, Time:=NULL] %>%
    .[, `:=`(donor=factor(donor),
             cell_type=factor(cell_type), 
             condition=factor(condition, levels=conds), 
             time=factor(time), 
             expt=factor(expt),
             sample_status=factor(sample_status))]
  
temp <- merge(dat.386, dat.490, all=T) %>%
  as.data.table() 
dt <- merge(temp, dat.537, all=T) %>%
  as.data.table() %>%
  .[, `:=`(donor=factor(donor, levels=donors),
           cell_type=factor(cell_type), 
           condition=factor(condition, levels=conds), 
           time=factor(time), 
           expt=factor(expt, levels=expts),
           sample_status=factor(sample_status, levels=sample.statuses))]

return(dt)
}

###### MAIN ######

dat.386 <- readFiles(p=paste0(fcs.path, "exp_00386/")) %>%
  combineFiles386(fd=paste0(table.path, "ARI_HC_EXP_00386_Bcell_process_data_2024.csv")) %>%
  asinTransform() %>%
  scaleData()
dat.490 <- readFiles(p=paste0(fcs.path, "exp_00490/")) %>%
  combineFiles490a537(fd=paste0(table.path, "ARI_HC_EXP_00490_Bcell_process_data_2024.csv")) %>%
  asinTransform() %>%
  scaleData() 
dat.537 <- readFiles(p=paste0(fcs.path, "exp_00537/")) %>%
  combineFiles490a537(fd=paste0(table.path, "ARI_HC_EXP_00537_Bcell_process_data_2024.csv")) %>%
  asinTransform() %>%
  scaleData() 
dat <- exptDataMerge()
dat[, .N, by=donor]

#save processed single-cell data file
write_csv(dat, paste0(table.path, "bcell_processed_flow_data.csv"))

#optional check of data file
#check.dat <- fread(paste0(table.path, "bcell_processed_flow_data.csv"))

