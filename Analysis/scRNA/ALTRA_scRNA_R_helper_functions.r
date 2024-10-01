# define the centralized functions to run analysis in ALTRA project

######## miscellaneous ########
# This code defines a color palette for RA converters.
# It reads a CSV file containing subject colors and assigns them to the corresponding subject GUIDs.
# The function `scale_color_subj` is defined to create a manual color scale using the subject colors.

subj_colors_list <- read_csv(
    file.path(
        "/home/jupyter/github/ra-longitudinal/metadata",
        "ALTRA_longtidunal_subject_colors.csv"
    ),
    show_col_types = FALSE
)
subject_colors <- c(subj_colors_list$subj_colors)
names(subject_colors) <- subj_colors_list$subject.subjectGuid

scale_color_subj <- function(colors = subject_colors, ...) {
    ggplot2:::manual_scale(
        "color",
        values = colors,
        ...
    )
}

######## DA analysis  ########

# define the centralized functions to run analysis in ALTRA project

# DA analysis

#' Calculate Pseudo Frequencies and Add Metadata
#'
#' This function calculates pseudo frequencies by adding +1 to the counts of each cell type
#' within each sample, then joins additional metadata and calculates live frequencies based
#' on the adjusted counts. This approach helps in normalizing data, especially in cases where
#' some cell types might not be present in a sample, thus avoiding division by zero errors
#' and ensuring a more robust analysis.
#'
#' @param so_labels_fl Data frame containing cell type labels and sample IDs.
#' @param celltype_column Name of the column containing cell type labels.
#' @param sample_id_column Name of the column containing sample IDs.
#' @return A data frame with pseudo frequencies and additional metadata.
#' @examples
#' # Assuming `so_labels_fl` is your data frame, `cellType` is the cell type column,
#' # and `sampleID` is the sample ID column:
#' result <- calculate_pseudo_frequencies(so_labels_fl, "cellType", "sampleID")
#' @export
calculate_pseudo_frequencies <- function(so_labels_fl, celltype_column, sample_id_column, metadata_columns) {
    # Convert the input data frame to a tibble and group by sample ID
    so_labels_freq <- so_labels_fl %>%
        as_tibble() %>%
        group_by(!!sym(sample_id_column)) %>%
        # Calculate the total counts for each sample
        mutate(total_counts = n()) %>%
        ungroup() %>%
        # Convert cell type and sample ID columns to factors
        mutate(
            !!celltype_column := factor(!!sym(celltype_column)),
            !!sample_id_column := factor(!!sym(sample_id_column))
        ) %>%
        # Group by cell type and sample ID, then calculate counts for each group
        group_by(!!sym(celltype_column), !!sym(sample_id_column), .drop = FALSE) %>%
        summarise(counts = n(), .groups = "drop") %>%
        ungroup() %>%
        # Add 1 to each count to calculate pseudo counts
        mutate(pseudo_counts = counts + 1)

    # Calculate the total pseudo counts for each sample
    pseudo_total_counts <- so_labels_freq %>%
        group_by(!!sym(sample_id_column)) %>%
        summarise(pseudo_total_counts = sum(pseudo_counts), .groups = "drop")

    # Join the total pseudo counts back to the main data frame
    so_labels_freq <- so_labels_freq %>%
        left_join(pseudo_total_counts, by = sample_id_column)

    # Select the specified metadata columns and total counts, ensuring each sample ID is unique
    meta_keep <- so_labels_fl %>%
        as_tibble() %>%
        select(!!sym(sample_id_column), all_of(metadata_columns)) %>%
        distinct(!!sym(sample_id_column), .keep_all = TRUE)

    # Join the metadata back to the main data frame
    so_labels_freq_meta <- so_labels_freq %>%
        left_join(meta_keep, by = sample_id_column) %>%
        # Calculate the live frequency for each cell type within each sample
        group_by(!!sym(sample_id_column)) %>%
        mutate(frequency = counts / sum(counts))

    # Calculate the pseudo frequency live for each cell type within each sample
    so_labels_freq_meta <- so_labels_freq_meta %>%
        group_by(!!sym(sample_id_column)) %>%
        mutate(pseudo_frequency = pseudo_counts / sum(pseudo_counts)) %>%
        ungroup()
    return(so_labels_freq_meta)
}


#' Perform CLR Transformation on Frequency Data
#'
#' This function applies a centered log-ratio (CLR) transformation to frequency data
#' from a given frequency table. It's designed to work with compositional data, making
#' the data more suitable for statistical analysis by addressing the issue of compositional
#' closure. The function allows specifying the columns for frequencies, cell types, and
#' sample IDs, and it returns a table with the original data and the CLR-transformed frequencies.
#'
#' @param freq_table A data frame containing the frequency data along with cell type labels
#' and sample IDs.
#' @param freq_col The name of the column in `freq_table` that contains the frequency data.
#' Defaults to 'frequency_live'.
#' @param celltype_col The name of the column in `freq_table` that contains the cell type labels.
#' Defaults to 'labels'.
#' @param sample_col The name of the column in `freq_table` that contains the sample IDs.
#' Defaults to 'sample_id'.
#' @param clr_name The name to be given to the column containing the CLR-transformed frequencies.
#' Defaults to 'clr'.
#' @return A data frame similar to `freq_table` but with an additional column for the CLR-transformed
#' frequencies.
#' @examples
#' # Assuming `freq_data` is your frequency table:
#' clr_transformed <- FreqClrTran(freq_data, "frequency_live", "labels", "sample_id", "clr")
#' @export
FreqClrTran <- function(freq_table, freq_col = "frequency_live", celltype_col = "labels", sample_col = "sample_id",
                        clr_name = "clr") {
    # Load the 'compositions' package for CLR transformation
    require("compositions")

    # Select and reshape the frequency data for CLR transformation
    freq <- freq_table %>%
        select(all_of(c(freq_col, celltype_col, sample_col))) %>%
        pivot_wider(
            id_cols = all_of(sample_col), names_from = all_of(celltype_col),
            values_from = all_of(freq_col), values_fill = 0
        )

    # Convert the reshaped data to a matrix format suitable for CLR transformation
    freq_mx <- freq %>%
        select(-all_of(sample_col)) %>%
        as.matrix()
    rownames(freq_mx) <- freq[[sample_col]]

    # Perform CLR transformation using the 'compositions' package
    freq_clr <- compositions::clr(freq_mx) %>%
        as_tibble(rownames = sample_col) %>%
        pivot_longer(cols = c(-all_of(sample_col)), names_to = celltype_col, values_to = clr_name)

    # Join the CLR-transformed frequencies back with the original frequency table
    freq_meta_clr <- freq_table %>% full_join(freq_clr, by = c(sample_col, celltype_col))

    # Return the combined data frame with CLR-transformed frequencies
    return(freq_meta_clr)
}

# Function: RunGlmFreq
# Description: Runs a generalized linear model (GLM) for frequency (percentage) data.
# Parameters:
#   - freq_table: A data frame containing the frequency table data.
#   - formula: A character string specifying the formula for the GLM.
#   - cell_type: A character string specifying the cell type.
# Returns:
#   - A data frame containing the tidy results of the GLM with an additional column for the cell type.

RunGlmFreq <- function(freq_table, formula, cell_type) {
    glm_res <- broom::tidy(stats::glm(as.formula(formula), data = freq_table)) %>%
        mutate(cell_type)
    return(glm_res)
}

# Function: RunLmeFreq
# Description: Runs a linear mixed model for time to conversion.
# Parameters:
#   - freq_table: A data frame containing the frequency table data.
#   - formula: A character string specifying the formula for the linear mixed model.
#   - celltype: A character string specifying the cell type being analyzed.
# Returns: A data frame containing the results of the linear mixed model analysis.

RunLmeFreq <- function(freq_table, formula, celltype) {
    message(paste("running", celltype))
    lme_res <- tidy(lmerTest::lmer(as.formula(formula), data = freq_table)) %>%
        mutate(celltype = celltype)
    return(lme_res)
}


######## DEG analysis  ########
# Function to create pseudo SingleCellExperiment objects for each cell type in AIM1
# Arguments:
#   - countfile_table: A data frame containing information about count files and metadata files for each cell type
# Returns:
#   - A list of pseudo SingleCellExperiment objects, with each object representing a cell type
makePseudoSE <- function(countfile_table) {
    require("SummarizedExperiment")
    pb_se <- lapply(1:nrow(countfile_table), function(x) {
        # load data for one cell type
        pb_counts <- fread(file.path(data_path, countfile_table$counts_file[x])) %>% rename("V1" = "index")
        # load the metadata
        pb_meta <- fread(file.path(data_path, countfile_table$meta_file[x])) %>% rename("V1" = "index")
        # check the index
        stopifnot(all(pb_counts$index == pb_meta$index))
        # make count matrix
        pb_counts_mx <- pb_counts %>%
            select(-index) %>%
            as.data.frame() %>%
            t()
        # make summarized experiment
        pb_se <- SummarizedExperiment(assays = list(counts = pb_counts_mx), colData = pb_meta)
    })
    names(pb_se) <- countfile_table$cell_type
    return(pb_se)
}

# RunDeseq function performs DESeq2 analysis on a list of SingleCellExperiment objects.
#
# Args:
#   se_list: A list of SingleCellExperiment objects.
#   formula: A formula specifying the design for DESeq2 analysis.
#
# Returns:
#   A list of DESeqDataSet objects, each representing the DESeq2 analysis result for a SingleCellExperiment object in se_list.
RunDeseq_list <- function(se_list, formula) {
    deseq_list <- lapply(seq_along(se_list), function(i) {
        require("DESeq2")
        pb <- se_list[[i]]
        cell_type <- names(se_list)[i]
        message(paste("run Deseq2 for cell type", cell_type))
        colnames(pb) <- colData(pb)$sample.sampleKitGuid
        pb_des <- DESeqDataSet(pb, design = as.formula(formula))
        pb_des <- DESeq(pb_des)
        return(pb_des)
    })
    names(deseq_list) <- names(se_list)
    return(deseq_list)
}

#' Add normalized counts to the SingleCellExperiment object
#'
#' This function adds normalized counts to a SingleCellExperiment object using the provided method.
#'
#' @param pseudo A SingleCellExperiment object to which normalized counts will be added.
#' @param dds A DESeqDataSet object containing the raw counts.
#' @param method The method to use for normalization. Can be one of "vst", "rlog", or "pseudocounts".
#' @param assay_name The name of the assay in the SingleCellExperiment object to store the normalized counts.
#'
#' @return The modified SingleCellExperiment object with normalized counts added.
#'
#' @examples
#' # Create a SingleCellExperiment object
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#'
#' # Add normalized counts using the "vst" method
#' sce <- addNormcounts(sce, dds, method = "vst")
#'
#' # Add normalized counts using the "rlog" method
#' sce <- addNormcounts(sce, dds, method = "rlog")
#'
#' # Add normalized counts using the "pseudocounts" method
#' sce <- addNormcounts(sce, dds, method = "pseudocounts")
#' @export
addNormcounts <- function(pseudo, dds, method = c("vst"),
                          assay_name = "normalized_counts") {
    estimateSizeFactors(dds)
    require("DESeq2")
    if (method == "vst") {
        norm_dds <- varianceStabilizingTransformation(dds, blind = FALSE)
    } else if (method == "rlog") {
        norm_dds <- rlog(dds, blind = FALSE)
    } else if (method == "pseudocounts") {
        norm_dds <- normTransform(dds)
    } else {
        stop("method should be one of vst/rlog/pseudocounts.")
    }
    norm_counts <- assay(norm_dds)
    colnames(norm_counts) <- NULL
    assay(pseudo, assay_name, withDimnames = FALSE) <- norm_counts
    return(pseudo)
}

# addNormcountsList function
#
# This function adds normalized counts from DESeq2 to a pseudobulk object.
#
# Parameters:
# - se_list: A list of pseudobulk objects.
# - deseq_list: A list of DESeq2 objects.
# - method: The method used for normalization (default is "vst").
# - assay_name: The name of the assay to store the normalized counts (default is "normalized_counts").
#
# Returns:
# A list of pseudobulk objects with added normalized counts.
#
# Example usage:
# se_list <- list(pseudobulk1, pseudobulk2, pseudobulk3)
# deseq_list <- list(deseq1, deseq2, deseq3)
# result <- addNormcountsList(se_list, deseq_list, method = "vst", assay_name = "normalized_counts")
#
addNormcountsList <- function(se_list, deseq_list, method = "vst",
                              assay_name = "normalized_counts") {
    se_norm_list <- lapply(names(se_list), function(i) {
        se <- se_list[[i]]
        cell_type <- i
        message(paste("Add normalized counts for", cell_type))
        se_list[[i]] <- addNormcounts(
            pseudo = se_list[[i]],
            dds = deseq_list[[i]],
            method = method, assay_name = assay_name
        )
        return(se_list[[i]])
    })
    names(se_norm_list) <- names(se_list)
    return(se_norm_list)
}

######## DEG analysis for linear mixed modeling ########
#' DEseqNormCount Function
#' This function takes a list of SingleCellExperiment objects and performs count normalization using the DESeq2 package.
#' It applies variance stabilization transformation (VST) to the counts, which helps to stabilize the variance and make the data more suitable for downstream analysis.
#' The function handles different biological sexes separately, as they may have major differences in the dataset.
#
#' @param se_list A list of SingleCellExperiment objects, each representing a single-cell RNA-seq dataset for a specific cell type.
#' @param cores The number of cores to be used for parallel processing. Default is 1.
#' @return A list of normalized SingleCellExperiment objects.
#' @import DESeq2
#' @import parallel
#' @export
DEseqNormCount <- function(se_list, cores = 1) {
    require("DESeq2")
    require("parallel")
    deseq_list <- mclapply(1:length(se_list), function(i) {
        se <- se_list[[i]]
        cell_type <- names(se_list)[i]
        message(paste("Processing", cell_type))
        colnames(se) <- colData(se)$sample.sampleKitGuid
        colData(se)$subject.biologicalSex <- factor(colData(se)$subject.biologicalSex)
        # # do variance stablization in sex differently as different sex have major difference in the dataset
        # if (length(unique(colData(se)$subject.biologicalSex)) >= 2) {
        #     pb_des <- DESeqDataSet(se, design = as.formula("~ subject.biologicalSex"))
        #     vsd <- varianceStabilizingTransformation(pb_des, blind = FALSE)
        # } else {
        #     pb_des <- DESeqDataSet(se, design = as.formula("~ 1"))
        #     vsd <- varianceStabilizingTransformation(pb_des, blind = TRUE)
        # }
        pb_des <- DESeqDataSet(se, design = as.formula("~ 1"))
        vsd <- varianceStabilizingTransformation(pb_des, blind = TRUE)
        assay(se, "vst") <- assay(vsd)
        return(se)
    }, mc.cores = cores)
    names(deseq_list) <- names(se_list)
    return(deseq_list)
}


# # run linear mixed model in lme4 directly
# # test in one gene
#' Run a linear mixed-effects model for one gene
#'
#' This function runs a linear mixed-effects model for a given gene using the lmerTest package.
#'
#' @param gene The gene of interest.
#' @param pb_des The data object containing the count data and other variables.
#' @param formula The formula specifying the model.
#' @param assays The type of assay to use for the count data (default is "vst").
#'
#' @return The fitted linear mixed-effects model.
#'
#' @examples
#' # Run the model for gene "MY_GENE"
#' model <- lme4_model_gene_norm(gene = "MY_GENE", pb_des = my_data, formula = y ~ x1 + x2)
#'
#' @import lmerTest
#' @importFrom dplyr select mutate
#' @importFrom SingleCellExperiment colData assay
#' @export
lme4_model_gene_norm <- function(gene, pb_des, formula, celltype_col, assays = "vst") {
    # extract count data
    count_data <- colData(pb_des) %>%
        as.data.frame() %>%
        select(
            all_of(c(
                celltype_col, "subject.subjectGuid", "sample.sampleKitGuid", "age_conv",
                "bmi_conv", "file.batchID", "Status_Long",
                "subject.biologicalSex", "days_to_conversion", "psbulk_n_cells"
            ))
        ) %>%
        mutate(yrs_to_conversion = days_to_conversion / 365)
    count_data[, "count"] <- assay(pb_des, assays)[gene, ]

    fit <- lmerTest::lmer(formula = formula, data = count_data)
    return(fit)
}


#' Run lme4 model for cell type
#'
#' This function runs the lme4 model using the DESeq2 normalized counts for a specific cell type.
#'
#' @param pb_des The DESeqDataSet object containing the normalized counts.
#' @param formula The formula specifying the model to be fitted.
#' @param assays The assays to be used for the model (default is "vst").
#' @param celltype_col The column in the colData specifying the cell type (default is "pred_manual").
#' @param cores The number of cores to be used for parallel processing (default is 62).
#' @param output_path The path where the output files will be saved.
#' @param proj_name The name of the project.
#'
#' @return A data frame containing the results of the lme4 model for each gene.
#'
#' @examples
#' runLme4Celltype(pb_des, formula, assays = "vst", celltype_col = "pred_manual", cores = 62, output_path = output_path, proj_name = proj_name)
runLme4Celltype <- function(pb_des, formula, assays = "vst", celltype_col = "pred_manual", cores = 62, output_path = output_path, proj_name = proj_name) {
    # load data
    cell_type <- colData(pb_des)[, celltype_col] %>% unique()
    message(paste("Run lme4 model deseq2 norm counts for cell type", cell_type))
    # load data
    message(paste("FIt glmm model  deseq2 norm counts for ", as.character(length(rownames(pb_des))), "genes"))
    # run a loop for all genes
    fit_lme4_norm <- mclapply(rownames(pb_des), lme4_model_gene_norm,
        pb = pb_des, assays = assays,
        formula = formula, mc.cores = cores, celltype_col = celltype_col
    )
    names(fit_lme4_norm) <- rownames(pb_des)
    # save model object
    saveRDS(fit_lme4_norm, file.path(output_path, paste0(
        proj_name, "_", cell_type,
        "_lme4_", assays, "_models.rds"
    )))
    # extract the results
    fit_lme4_res <- lapply(names(fit_lme4_norm), function(x) {
        res <- fit_lme4_norm[[x]] %>%
            tidy() %>%
            mutate(gene = x)
        return(res)
    }) %>%
        data.table::rbindlist() %>%
        mutate(celltype = cell_type)
    # save model object
    write_csv(fit_lme4_res, file.path(output_path, paste0(
        proj_name, "_", cell_type,
        "_lme4_", assays, "_models_results.csv"
    )))

    return(fit_lme4_res)
}


# This function calculates the q-values for a given data frame based on the p-values in a specified column.
# Args:
#   x: A data frame containing the p-values.
#   p_col: The name of the column in the data frame that contains the p-values. Default is 'p.value'.
# Returns:
#   The input data frame with an additional column 'q_values' containing the calculated q-values.
calQvalue <- function(x, p_col = "p.value") {
    require("qvalue")
    pvalues <- x %>% pull(.data[[p_col]])
    qvalues <- qvalue::qvalue(pvalues)$qvalues
    x <- x %>% mutate("q_values" = qvalues)
    return(x)
}


#' Run Gene Set Enrichment Analysis (GSEA) per cell type
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) on a list of genes for each cell type.
#'
#' @param logfc_list A data frame containing the gene expression data with columns for cell type, log fold change, gene, and other optional columns.
#' @param gmx A gene set matrix (GMT) file or a pre-loaded gene set matrix. If not provided, the function will load and initialize the pathway database using a default GMT file.
#' @param ct.col The column name in the logfc_list data frame that represents the cell type.
#' @param rank.col The column name in the logfc_list data frame that represents the log fold change.
#' @param gene.col The column name in the logfc_list data frame that represents the gene.
#' @param collapsePathways Logical value indicating whether to collapse pathways with similar gene sets into a single pathway.
#' @param ncores The number of cores to use for parallelization. If not provided, the function will use all available cores minus 3.
#'
#' @return A data frame containing the results of the GSEA analysis, including pathway information, p-values, adjusted p-values, normalized enrichment scores (NES), and leading edge genes.
#'
#' @import fgsea
#' @import dplyr
#' @import rlist
#' @import BiocParallel
#' @import data.table
#'
#' @examples
#' # Example usage of RunGSEACelltype function
#' logfc_list <- data.frame(
#'     cell_type = c("A", "A", "B", "B"),
#'     logfc = c(1.2, 0.8, 1.5, 0.5),
#'     gene = c("gene1", "gene2", "gene3", "gene4")
#' )
#' result <- RunGSEACelltype(logfc_list)
#' print(result)
#'
#' @export

RunGSEACelltype <- function(logfc_list, gmx = NULL, ct.col = "cell_type", rank.col = "logfc", gene.col = "gene",
                            collapsePathways = FALSE,
                            ncores = 1) {
    require(fgsea)
    # if no provided, Load and initialize pathway database
    if (is.null(gmx)) {
        gmxFile <- "/home/jupyter/data/reference/c2.cp.v7.2.symbols.gmt"
        colNames <- max(count.fields(file = gmxFile, sep = "\t"))
        colNames <- seq(from = 1, to = colNames)
        colNames <- as.character(colNames)
        gmx <- read.table(
            file = gmxFile,
            sep = "\t",
            quote = "\"",
            fill = TRUE,
            col.names = colNames,
            row.names = 1
        )
        gmx <- gmx[, -1]
        gmx <- apply(gmx, MARGIN = 1, FUN = function(x) {
            return(value = setdiff(unname(x), ""))
        })
        names(gmx) <- toupper(names(gmx))
    }

    # setup parallelization parameters
    if (is.null(ncores)) {
        ncores <- parallel::detectCores() - 3
    } else {
        (ncores <- ncores)
    }
    param <- BiocParallel::MulticoreParam(workers = ncores, progressbar = TRUE)

    # RUN GSEA per celltype
    celltypes <- unique(logfc_list %>% pull(.data[[ct.col]]))

    pLS <- lapply(celltypes, function(ct) {
        message(paste("run GSEA in", ct))

        # create rank list based on lowest to higest gene fold-change
        rnkDF <- logfc_list %>%
            dplyr::filter(.data[[ct.col]] == ct) %>%
            filter(!is.na(.data[[rank.col]])) %>%
            dplyr::arrange(.data[[rank.col]])
        rnk <- rnkDF %>%
            pull(.data[[rank.col]]) %>%
            as.numeric()
        names(rnk) <- rnkDF %>% pull(.data[[gene.col]])
        message(paste("run GSEA in", length(rnk), "genes"))

        # run GSEA by parallelization
        fgseaRes <- fgsea::fgsea(
            pathways = gmx,
            stats = rnk,
            minSize = 10,
            maxSize = 500, nPermSimple = 10000,
            BPPARAM = param
        )

        # filter on pathways <0.05 adjusted p-value
        fgseaRes_tb <- fgseaRes %>%
            as.data.frame() %>%
            dplyr::filter(padj < 0.05) %>%
            dplyr::select(pathway, pval, padj, NES, leadingEdge) %>%
            dplyr::arrange(desc(NES)) %>%
            dplyr::mutate(celltype = ct)
        # if only keep the main pathway
        if (collapsePathways) {
            collapsedPathways <- fgsea::collapsePathways(
                fgseaRes[order(pval)][padj < 0.05],
                gmx, rnk
            )
            mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                order(-NES), pathway
            ]
            fgseaRes_tb <- fgseaRes_tb %>% dplyr::filter(pathway %in% mainPathways)
        }

        return(value = fgseaRes_tb)
    })

    pathwayDF <- data.table::rbindlist(pLS)
    pathwayDF$leadingEdge <- vapply(pathwayDF$leadingEdge,
        paste,
        collapse = ", ",
        character(1L)
    )

    # make plotting data frame
    plotDF <- pathwayDF %>%
        mutate(
            group = ifelse(NES > 0, "up", "down"),
            pID = c(1:length(pathway))
        ) %>%
        group_by(pID) %>%
        mutate(lesize = length(unlist(strsplit(leadingEdge, ",")))) %>%
        as_tibble()


    # determine pathway size
    gsSize <- data.frame(gsize = sapply(gmx, function(x) length(x))) %>%
        rownames_to_column(var = "pathway")

    # calculate propotion of genes enriched (#leading edge genes/size of pathway)
    plotDF <- plotDF %>%
        mutate(
            gsize = gsSize$gsize[match(pathway,
                table = gsSize$pathway
            )],
            propGenes = (lesize / gsize) * 100
        )

    return(plotDF)
}


#' Normalize counts using DESeq and variance stabilizing transformation
#'
#' This function takes a list of SummarizedExperiment objects, applies DESeq normalization
#' based on the biological sex of the subjects, and then applies a variance stabilizing
#' transformation. If the biological sex is not a variable (less than 2 unique values),
#' it falls back to a basic normalization without considering biological sex.
#'
#' @param aim3_pb_fl A list of SummarizedExperiment objects.
#' @return A list of SummarizedExperiment objects with an added assay named 'vst' containing
#'         the variance stabilized transformed values.
#' @examples
#' # Assuming `aim3_pb_fl` is your list of SummarizedExperiment objects:
#' normalized_data <- normalize_counts(aim3_pb_fl)
#' @export
normalize_counts <- function(aim3_pb_fl) {
    aim3_pb_deseq <- lapply(aim3_pb_fl, function(x) {
        colnames(x) <- colData(x)$sample.sampleKitGuid
        colData(x)$subject.biologicalSex <- factor(colData(x)$subject.biologicalSex)
        if (length(unique(colData(x)$subject.biologicalSex)) >= 2) {
            pb_des <- DESeqDataSet(x, design = as.formula("~ subject.biologicalSex"))
            vsd <- varianceStabilizingTransformation(pb_des, blind = FALSE)
        } else {
            pb_des <- DESeqDataSet(x, design = as.formula("~ 1"))
            vsd <- varianceStabilizingTransformation(pb_des, blind = TRUE)
        }
        assay(x, "vst") <- assay(vsd)
        return(x)
    })
    return(aim3_pb_deseq)
}


#' Get expression data for differentially expressed genes in single-cell RNA-seq data
#'
#' This function takes a list of differentially expressed genes (DEGs) and a list of single-cell experiments (SCEs)
#' and returns a data frame containing the expression data for each DEG in each cell type.
#'
#' @param deg_list A data frame containing the DEGs, with columns "cell_type" and "gene".
#' @param SE_list A list of single-cell experiments, where each element corresponds to a different cell type.
#' @param assay The assay type to use for extracting gene expression data (default is "normalized_counts").
#'
#' @return A data frame containing the expression data for each DEG in each cell type, with columns "cell_type",
#' "gene", and "normalized_counts".
#'
#' @examples
#' # Example usage
#' deg_list <- data.frame(cell_type = c("T cell", "B cell"), gene = c("CD3D", "CD19"))
#' SE_list <- list(T_cell = pb_sce_T_cell, B_cell = pb_sce_B_cell)
#' gene_exprs <- GetExprsSE(deg_list, SE_list)
#'
#' @import data.table
#' @import dplyr
#' @import SingleCellExperiment
GetExprsSE <- function(deg_list, SE_list, assay = "normalized_counts") {
    gene_exprs <- lapply(1:nrow(deg_list), function(i) {
        cell_type_test <- deg_list$cell_type[i]
        gene_test <- deg_list$gene[i]
        pb <- SE_list[[cell_type_test]]
        # Convert the data into sce object for plotting
        pb_sce <- as(pb, "SingleCellExperiment")
        # Get expression matrix
        gex_exprs <- assay(pb_sce, assay)[gene_test, ]
        # Set up the column annotation
        pb_df <- colData(pb_sce) %>%
            as_tibble() %>%
            mutate(
                cell_type = cell_type_test,
                gene = gene_test,
                normalized_counts = gex_exprs
            )
        return(pb_df)
    }) %>%
        data.table::rbindlist() %>%
        mutate(celltype_gene = paste0(cell_type, ": ", gene))
    return(gene_exprs)
}


#' PlotExprsLgPair Function
#'
#' This function generates a line plot and a violin plot for gene expression data.
#'
#' @param gex_lg A data frame containing gene expression data for the large cohort.
#' @param gex_conv A data frame containing gene expression data for the conversion cohort.
#' @param celltype_plot The cell type to plot.
#' @param gene_plot The gene to plot.
#' @param subject_colors A vector of colors for subject IDs.
#' @param group_colors A vector of colors for group IDs.
#' @param fig_path The path to save the generated plot.
#' @param proj_name The name of the project.
#'
#' @return A combined plot of the line plot and the violin plot.
#'
#' @import ggplot2
#' @import dplyr
#' @import ggpubr
#'
#' @examples
#' PlotExprsLgPair(gex_lg, gex_conv, "T cell", "CD3D",
#'     subject_colors = cluster_colors, group_colors = ari_colors,
#'     fig_path = "/path/to/save/plot", proj_name = "Project1"
#' )
PlotExprsLgPair <- function(gex_lg, gex_conv, celltype_plot, gene_plot,
                            subject_colors = cluster_colors, group_colors = ari_colors, fig_path, proj_name) {
    library(ggplot2)
    library(dplyr)
    library(ggpubr)

    p1 <- gex_lg %>%
        filter(cell_type == celltype_plot & gene == gene_plot) %>%
        ggplot(aes(x = days_to_conversion, y = normalized_counts)) +
        geom_line(aes(color = subject.subjectGuid)) +
        geom_point(aes(color = subject.subjectGuid, shape = Sex)) +
        stat_smooth(
            color = "#8F7700FF", fill = "#8F7700FF",
            method = "loess"
        ) +
        labs(title = paste(celltype_plot, gene_plot), x = "Days to Conversion", y = "Normalized expression") +
        scale_color_subj() +
        theme_minimal() +
        theme(
            text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 10, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            plot.title = element_text(size = 10, face = "bold"),
            legend.position = "top"
        ) +
        scale_shape_manual(values = c(1, 2)) +
        guides(color = FALSE)

    # only plot the paired samples
    gex_conv_fl <- gex_conv %>%
        filter(cell_type == celltype_plot & gene == gene_plot)
    pair_nums <- gex_conv_fl %>%
        group_by(subject.subjectGuid) %>%
        tally()
    pair_subs <- pair_nums %>%
        filter(n >= 2) %>%
        pull(subject.subjectGuid) %>%
        unique()
    p2 <- gex_conv_fl %>%
        filter(subject.subjectGuid %in% pair_subs) %>%
        ggplot(aes(x = status, y = normalized_counts)) +
        geom_violin(aes(), alpha = 0.6) +
        geom_point(aes(group = subject.subjectGuid, color = subject.subjectGuid, shape = Sex)) +
        geom_line(aes(group = subject.subjectGuid, color = subject.subjectGuid)) +
        labs(title = "", x = "", y = "") +
        scale_color_subj() +
        scale_fill_manual(values = group_colors[2:3]) +
        theme_minimal() +
        theme(
            text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 10, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            plot.title = element_text(size = 10, face = "bold"),
            legend.position = "top", axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
        ) +
        scale_shape_manual(values = c(1, 2)) +
        guides(color = FALSE, fill = FALSE)

    p_combine <- ggarrange(p1, p2, common.legend = TRUE)
    ggsave(file.path(fig_path, paste0(
        proj_name, celltype_plot,
        "_", gene_plot, "_expression.pdf"
    )), width = 6, height = 4)
    return(p_combine)
}
