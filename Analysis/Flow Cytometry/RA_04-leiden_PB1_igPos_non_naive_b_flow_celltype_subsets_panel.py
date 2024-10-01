#!/usr/bin/python3
 
import scyan as sy
import os
import glob
import anndata
import re
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import scanpy as sc
import scanpy.external as sce
import copy
from concurrent.futures import ProcessPoolExecutor

print(sy.__version__)
print(anndata.__version__)

#### ==== parameter set up ==== ####
# define the working path
panel = "PB1"
# define l1 cell types to subset
cell_types = ['total_b_cells_igPos_non_naive_b']

data_path='/home/jupyter/projects/pre-ra/flow/raw-data/' + panel + '/labelled-expr/cache/'
fig_path = '/home/jupyter/projects/pre-ra/flow/02-clustering/results/' + panel + '_igPos_non_naive_b_celltype_subsample'  + "/"
proj_name = 'pre-ra_flow_clustering_' + panel
output_path = '/home/jupyter/projects/pre-ra/flow/02-clustering/data/' +panel + '/'


if not os.path.exists(fig_path): os.makedirs(fig_path)
    
if not os.path.exists(output_path): os.makedirs(output_path)

# define scanpy verbose levels
sc.settings.verbosity = 3
sc.settings.figdir = fig_path
sc.settings.n_jobs = -1

#### ==== FUNCTIONS ==== ####
def read_anndata_files(file_tuples):
    """
    Read Anndata objects from H5AD files and store them in a dictionary with custom names.

    Parameters:
        file_tuples (list of tuples): List of tuples where each tuple contains filename and desired name.

    Returns:
        dict: Dictionary containing Anndata objects with custom names.
    """
    anndata_dict = {}
    for filename, name in file_tuples:
        anndata_obj = anndata.read_h5ad(filename)
        anndata_dict[name] = anndata_obj
    return anndata_dict

def save_anndata_list(anndata_list, filenames):
    """
    Save a list of Anndata objects to separate H5AD files.

    Parameters:
        anndata_list (list): List of Anndata objects to be saved.
        filenames (list): List of filenames to save Anndata objects.

    Returns:
        None
    """
    # zip() combines two lists into a single tuple
    for anndata_obj, filename in zip(adata_downsampled_celltypes, filenames):
        print(adata_downsampled_celltypes[anndata_obj])
        adata_downsampled_celltypes[anndata_obj].write_h5ad(filename)
        

# make a function to find files
def get_filepaths_with_glob(root_path: str, file_regex: str):
    return glob.glob(os.path.join(root_path, file_regex))

#### parallelized leiden clustering
def run_leiden(adata, resolution, key_added):
    # Make a copy of adata for thread safety
    adata_copy = copy.deepcopy(adata)
    adata_clustering = sc.tl.leiden(adata_copy, resolution=resolution, key_added=key_added, copy=True)
    return adata_clustering.obs
 
def run_leiden_parallel(adata, tasks):
    with ProcessPoolExecutor(max_workers=5) as executor:
        # Make deep copies of adata for each task to ensure thread safety
        futures = [executor.submit(run_leiden, copy.deepcopy(adata), resolution, key_added) for resolution, key_added in tasks]
        
        results = [future.result() for future in futures]
 
    # Assign the results back to the original AnnData object
    for result, (_, key_added) in zip(results, tasks):
        adata.obs[key_added] = result[key_added]
 
    return adata
#### ==== READ IN ==== ####

filenames = get_filepaths_with_glob(output_path, "adata_leiden_scaled_harmonized_umap_downsmpl_celltypes_" + cell_types[0]+ "*.h5ad")  
print(filenames)

# Extract the cell types between "celltypes_" and "PT1.h5ad"
cell_types = [path.split('celltypes_')[1].split('_' +panel+'.h5ad')[0] for path in filenames]
print(cell_types)

## turn filenames and cell types into tuple
file_tups = list(zip(filenames, cell_types))
print(file_tups)

adata_downsampled_celltypes = read_anndata_files(file_tups)
print(adata_downsampled_celltypes)

#### ==== LEIDEN ==== ####
''''
leiden_fig_path = '/home/jupyter/projects/pre-ra/flow/02-clustering/results/' + proj_name + "/" + panel + "_leiden/"
print(leiden_fig_path)
if not os.path.exists(leiden_fig_path): os.makedirs(leiden_fig_path)


for label, sub_adata in adata_downsampled_celltypes.items():
    print("Key:", label)
    print("Value:", sub_adata)
    
    tasks = [(2, "leiden_res_2"),(2.5, "leiden_res_2.5"),(3, "leiden_res_3")]
    leiden_res = ["leiden_res_2", "leiden_res_2.5",  "leiden_res_3"]

    sub_adata = run_leiden_parallel(sub_adata, tasks)

    print('leiden clustering completed...')
    
    p1 = sy.plot.umap(sub_adata, color= leiden_res ,return_fig = True, size = .5, ncols = 3,legend_loc="on data")
    p1.savefig(leiden_fig_path + label   + "_umap_leiden_resolutions_" + panel + ".png",  dpi=400, bbox_inches='tight')
    
    for res in leiden_res: 
        p1 = sy.plot.umap(sub_adata, color= res ,return_fig = True, size = .5)
        p1.savefig(leiden_fig_path + label   + "_umap_" + res + "_" +panel + ".png",dpi=400, bbox_inches='tight')
            
    # store back to dict
    adata_downsampled_celltypes[label] = sub_adata
    
#### ==== SAVE ==== ####
filenames = [output_path + "adata_leiden_scaled_harmonized_umap_downsmpl_celltypes_" + cell_type + "_" + panel + ".h5ad" for cell_type in cell_types]
print(filenames)

save_anndata_list(adata_downsampled_celltypes, filenames)

''''

#### ==== LEIDEN SPEED UP ==== ####
# create separate leiden folder
leiden_fig_path = fig_path +  "/" + panel + "_leiden/"
print(leiden_fig_path)
if not os.path.exists(leiden_fig_path): os.makedirs(leiden_fig_path)


keys_to_subset = cell_types

subset_data = {key: value for key, value in adata_downsampled_celltypes.items() if key in keys_to_subset}
print(subset_data)


for label, sub_adata in subset_data.items():
    print("Key:", label)
    print("Value:", sub_adata)
    
    tasks = [(1.0, "leiden_res_1.0")]
    leiden_res = [ "leiden_res_1.0"]

    sub_adata = run_leiden_parallel(sub_adata, tasks)

    print('leiden clustering completed...')
    
    p1 = sy.plot.umap(sub_adata, color= ['leiden_res_2', 'leiden_res_2.5', 'leiden_res_3', 'leiden_res_1.0'] ,return_fig = True, size = .5, ncols = 3,legend_loc="on data")
    p1.savefig(leiden_fig_path + label   + "_umap_leiden_resolutions_" + panel + ".png",  dpi=400, bbox_inches='tight')
    
    for res in leiden_res: 
        p1 = sy.plot.umap(sub_adata, color= res ,return_fig = True, size = .5,legend_loc="on data")
        p1.savefig(leiden_fig_path + label   + "_umap_" + res + "_" +panel + ".png",dpi=400, bbox_inches='tight')
            
     # store back to dict
    adata_downsampled_celltypes[label] = sub_adata
    
    sub_adata.write_h5ad(output_path + "adata_leiden_scaled_harmonized_umap_downsmpl_celltypes_" + label + "_" + panel + ".h5ad")
    
    
    