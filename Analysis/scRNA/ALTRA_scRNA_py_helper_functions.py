def generate_pseudobulk_data(adata_pb, output_path, cell_type_col, proj_name, group='Status_Xsec'):
    """
    Generate pseudobulk data for each cell type in the provided AnnData object.

    Parameters:
        adata_pb (AnnData): The AnnData object containing the single-cell RNA-seq data.
        output_path (str): The path to the directory where the pseudobulk data will be saved.
        cell_type_col (str): The name of the column in the `obs` attribute of `adata_pb` that contains the cell type information.
    Returns:
        None
    """
    import decoupler as dc
    import pandas as pd
    for cell_type in adata_pb.obs[cell_type_col].unique():
        print("generating pseudobulk data for " + cell_type)
        adata_pb_sub = adata_pb[adata_pb.obs[cell_type_col]
                                == cell_type].copy()
        sample_number = len(adata_pb_sub.obs['sample.sampleKitGuid'].unique())
        # # filter pseudobulk data by proportion expression
        adata_pb_sub = adata_pb_sub[:, dc.filter_by_prop(
            adata_pb_sub, min_prop=0.1, min_smpls=5)].copy()
        # filter pseudobulk data by gene counts
        # similiar parameters used by filter_by_expr implemented in edgeR
        genes_keep = dc.filter_by_expr(
            adata_pb_sub, group=group, min_count=10, min_total_count=15)
        print("remaining " + str(len(genes_keep)) + ' genes in ' + cell_type)
        adata_pb_sub = adata_pb_sub[:, genes_keep].copy()
        # output counts and metadata file
        # save the counts data
        adata_pb_sub_counts = pd.DataFrame(index=adata_pb_sub.obs.index, columns=adata_pb_sub.var.index,
                                           data=adata_pb_sub.layers['counts'])
        adata_pb_sub_counts.to_csv(output_path+proj_name + '_' +
                                   cell_type + '_psbulk_counts.tsv', sep="\t")
        # save the meta data
        adata_pb_sub.obs.to_csv(output_path + proj_name + '_' + cell_type
                                + '_psbulk_metadata.csv')
