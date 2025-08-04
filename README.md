# ALTRA-manuscript
This repository contains Jupyter notebooks and scripts for the ALTRA manuscript data analysis and figure generation.

## Notebooks Overview

### Root Directory

- **data_preprocessing.ipynb**  
    Cleans raw data, performs normalization, and applies quality control for downstream analysis.

- **statistical_analysis.ipynb**  
    Runs statistical tests and fits models to identify significant biological associations.

- **figure_generation.ipynb**  
    Creates main manuscript figures from processed data.

- **supplementary_analysis.ipynb**  
    Performs additional analyses, including sensitivity checks and alternative methods.

### `Analysis` Folder

- **Flow_Cytometry/RA_01-clustering_flow_panels_PB1.ipynb**
- **Flow_Cytometry/RA_01-clustering_flow_panels_PT1.ipynb**
- **Flow_Cytometry/RA_02-sample_selection_PB1_flow_celltype_subsets_panels.ipynb**
- *(other Flow_Cytometry notebooks...)*
- **scRNA/ALTRA_scRNA_py_helper_functions.py**
- *(other scRNA notebooks and scripts...)*

### `Figures` Folder

#### Figure1
- **Figure1/Fig_S3E.ipynb**
- **Figure1/Figure_1B.ipynb**
- *(Autoantibody levels across disease status)*
- **Figure1/Figure_1F.ipynb**
- *(Differentially Expressed Genes between CON1 and ARI)*
- **Figure1/Figure_S1B.ipynb**
- *(Timepoints of longitunidal RA Converters to clinically significant RA)*
- **Figure1/Figure_S1C.ipynb**
- *(RF levels across the groups)*
- **Figure1/Figure_S2E.ipynb**
- *(Distribution of disease status across the OLINK clusters)*
- **Figure1/Figure_S3_A_B_C.ipynb**
- *(Celltype labelled single cell level UMPAs)*
- **Figure1/Figure_S3D.ipynb**
- *(Frequency changes overtime between CON1 and ARI)*

#### Figure2
- **Figure2/Figure_2C.ipynb**
- *(Frequency changes overtime and corresponding DEG counts in ARI)*
- **Figure2/Figure_2D.ipynb**
- *(DEGs for Longitudinal RA Converters)*
- **Figure2/Figure_2F_S5G.ipynb**
- *(paired Pseudobulk DEG analysis results)*
- **Figure2/Figure_2I.ipynb**
- *(IL1B+ monocyte frequency pre and post conversion in RA converters)*
- **Figure2/Figure_S5_A_C.ipynb**
- *((i) Autoanibodies levels in RA Converters)*
- *((ii) Longitudinal time points for CON2 subjects)*
- **Figure2/Figure_S5_D_E.ipynb**
- *(Flow cytometery :T cell clustering )*
- **Figure2/Figure_S5_F.ipynb**
- *(Flow cytometery : Frequency changes in CM/EM CD4 Tcells and Naive CD4 T cells)*
- **Figure2/Figure_SH.ipynb**
- *(DEGs in Monocytes )*

#### Figure3
- **Figure3/Figure_3D.ipynb**
- *(Frequency of Isotype population in B cells across Leiden clusters)*
- **Figure3/Figure_3G.ipynb**
- *(Flow cytometery :Frequency changes for ARI in Naive B cells)*
- **Figure3/Figure_S4F.ipynb**
- *(Frequency of Isotype population in B cells)*
- **Figure3/Figure_S5.ipynb**
- *(Proportion of Isotypes in B cells)*
- **Figure3/Figure_S6_G.ipynb**
- *(Key marker gene expression across Memory and Effector B cells)*
- **Figure3/Figure_S6_J.ipynb**
- *(Flow cytometery :B cell clustering)*

#### Figure4
- **Figure4/Figure_4A.ipynb**
- *(DEGs for CM CD4 T cells)*
- **Figure4/Figure_4B.ipynb**
- *(Gene module score plot for CM CD4 T cells)*
- **Figure4/Figure_4D_S7E.ipynb**
- *(Frequnecy for Non-naive CD4 T cells)*
- **Figure4/Figure_4G.ipynb**
- *(DEGs for Non-naive CD4 T cells at single cell level)*
- **Figure4/Figure_4_C_E_F_H_S6_S7.ipynb**
- *(Non-Naive CD4 T cells Plots- NMF projection)*
- **Figure4/Figure_S8C.ipynb**
- *(Mean gene expression of DEGs in Non-naive CD4 Tcells)*

#### Figure5
- **Figure5/Figure_5ABCD.ipynb**
- *(...)*
- **Figure5/Figure_5BCDE.ipynb**
- *(...)*
- **Figure5/Figure_5F.ipynb**
- *(DEGs in CD4 Naive T cells)*
- **Figure5/Figure_5G.ipynb**
- *(Gene module score plot for CD4 Naive T cells)*
- **Figure5/Figure_S11_ABEF.ipynb**
- *(Frequency for Core naive CD4 Tcells based on key DEGs)*
- **Figure5/Figure_S11D.ipynb**
- *(DEGs CD8 Naive T cells)*
- **Figure5/Figure_S11E.ipynb**
- *(Gene module score plot for CD8 Naive T cells)*

#### Figure6
- *(Figure6 notebooks...)*

#### Figure7
- *(Figure7 notebooks...)*

#### FigureS10
- *(FigureS10 notebooks...)*

## Usage

Clone the repository and run the notebooks in order. See `environment.yml` for dependencies.

For questions, contact the repository maintainer.
