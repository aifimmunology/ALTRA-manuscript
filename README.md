# ALTRA-manuscript
This repository contains Jupyter notebooks and scripts for the ALTRA manuscript data analysis and figures

## Notebooks Overview

### `Analysis` Folder

- **Flow_Cytometry**
- *Data cleaning, pre-processing ,Object creation and downstream analysis for Flow Cytometry data*
- **Intercellular_Bcell_Analysis**
- *Initial B cell data processing Metaclustering, isotype labeling and cleaning of processed B cell cytometry data*
- **OLINK_analysis**
- *Data cleaning, pre-processing and clustering analysis for Proteomics/OLINK data*
- **TEA-Seq_analysis**
- *Data cleaning, pre-processing ,Object creation and downstream analysis including MOFA and ACTAC seq for TEA-Seq data*
- **scRNA**
- *Data cleaning, pre-processing ,Object creation and downstream analysis including  Differential Abundance/Frequency , NMF analysis, Spectra analysis,Tfh17 validation analysis for scRNA data*
- NMF Analysis set up : https://github.com/yyoshiaki/NMFprojection/ 

### `Figures` Folder

#### Figure1
- **Figure1/Fig_S3E.ipynb**
- *Frequency changes in CM CD4 T cells and Core CD14 Monocytes between CON1/HC1 and ARI*
- **Figure1/Figure_1B.ipynb**
- *Autoantibody levels across disease status*
- **Figure1/Figure_1F.ipynb**
- *Differentially Expressed Genes between CON1/HC1 and ARI*
- **Figure1/Figure_S1B.ipynb**
- *Timepoints of longitunidal RA Converters to clinically significant RA*
- **Figure1/Figure_S1C.ipynb**
- *RF levels across the groups*
- **Figure1/Figure_S2E.ipynb**
- *Distribution of disease status across the OLINK clusters*
- **Figure1/Figure_S3_A_B_C.ipynb**
- *Celltype labelled single cell level UMPAs*
- **Figure1/Figure_S3D.ipynb**
- *Frequency changes overtime between CON1/HC1 and ARI*

#### Figure2
- **Figure2/Figure_2B.ipynb**
- *Intra donor variablity between CON2/HC2 and ARI*
- **Figure2/Figure_2C.ipynb**
- *Frequency changes overtime and corresponding DEG counts in ARI*
- **Figure2/Figure_2D.ipynb**
- *DEGs for Longitudinal RA Converters*
- **Figure2/Figure_2F_S5G.ipynb**
- *paired Pseudobulk DEG analysis results*
- **Figure2/Figure_2I.ipynb**
- *IL1B+ monocyte frequency pre and post conversion in RA converters*
- **Figure2/Figure_S4ABE.ipynb**
- *Gene Analysis reuslts using spectra model*
- **Figure2/Figure_S4CD.ipynb**
- *DEG analysis between RA converters , ARI and Healthy controls at Baseline*
- **Figure2/Figure_S5B.ipynb**
- *Clinical lab data analysis results of ALTRA cohort*
- **Figure2/Figure_S5_A_C.ipynb**
- *(i) Autoanibodies levels in RA Converters*
- *(ii) Longitudinal time points for CON2/HC2 subjects*
- **Figure2/Figure_S5_D_E.ipynb**
- *Flow cytometery :T cell clustering*
- **Figure2/Figure_S5_F.ipynb**
- *Flow cytometery : Frequency changes in CM/EM CD4 Tcells and Naive CD4 T cells*
- **Figure2/Figure_S5_H.ipynb**
- *DEGs in Monocytes*

#### Figure3
- **Figure3/Figure_3B.ipynb**
- *DEGs in Effector B cells*
- **Figure3/Figure_3C.ipynb**
- *Proportion of Isotypes in Effector B cells*
- **Figure3/Figure_3F.ipynb**
- *Flow Cytometry: Frequency changes in longitunal RA converters in Naive B cells*
- **Figure3/Figure_3I.ipynb**
- *Key marker gene expression between HC2 and ARI*
- **Figure3/Figure_S6A.ipynb**
- *Frequency changes in longitunal RA converters in Naive B cells*
- **Figure3/Figure_S6C.ipynb**
- *Key Isotype expression in Effector B cells*
- **Figure3/Figure_S6E.ipynb**
- *Frequency changes in HC2 in Effector B cells*
- **Figure3/Figure_S4_G.ipynb**
- *Expession of Isotypes across B cells*
- **Figure3/Figure_S6_H.ipynb**
- *Key marker gene expression across B cells*
- **Figure3/Figure_S6K.ipynb**
- *Flow cytometry: B cell Clustering*
#### Figure3_IC
- **Figure_3_IC/IC_flow_bcell_data_analysis.R**
- *Intracellular flow cytometry: Percentage cytokine-positive cells among naive B cells of ARI and HC2/CON2*
  
#### Figure4
- **Figure4/Figure_4A.ipynb**
- *DEGs for CM CD4 T cells*
- **Figure4/Figure_4B.ipynb**
- *Gene module score plot for CM CD4 T cells*
- **Figure4/Figure_4C_S7E.ipynb**
- *Frequnecy for Non-naive CD4 T cells*
- **Figure4/Figure_4E.ipynb**
- *Tfh/Tph gene score across clusters*
- **Figure4/Figure_4H.ipynb**
- *DEGs for Non-naive CD4 T cells at single cell level*
- **Figure4/Figure_4_DFGH_S6_S7.ipynb**
- *Non-Naive CD4 T cells Plots- NMF projection*
- **Figure4/Figure_S8CDEFG.ipynb**
- *Tfh17 signature in TEA-Seq data*
- **Figure4/Figure_S8H.ipynb**
- *Mean gene expression of DEGs in Non-naive CD4 T cells*
- **Figure4/Figure_S9BCD.ipynb**
- *Flow Cytometry : tfh17 clustering*

#### Figure5
- **Figure5/Figure_5ABCD.ipynb**
- *TEA-Seq: MOFA analysis related figures for CD4 Naive T cells*
- **Figure5/Figure_5BCDE.ipynb**
- *TEA-Seq: MOFA model results for CD4 Naive T cells*
- **Figure5/Figure_5F.ipynb**
- *DEGs in CD4 Naive T cells*
- **Figure5/Figure_5G.ipynb**
- *Gene module score plot for CD4 Naive T cells*
- **Figure5/Figure_S11_EF.ipynb**
- *Frequency for Core naive CD4 Tcells based on key DEGs*
- **Figure5/Figure_S11D.ipynb**
- *DEGs CD8 Naive T cells*
- **Figure5/Figure_S11E.ipynb**
- *Gene module score plot for CD8 Naive T cells*

#### Figure6_S12
- **Figure6_S12/Figure_6ADE_S12C.ipynb**
- *TEA-Seq : ATAC analysis figures for CD4 T cells*
- **Figure6_S12/Figure_6B_S12A.ipynb**
- *TEA-Seq : Frequency of CD4 T cells Clusters*
- **Figure6_S12/Figure_6C_S12B.ipynb**
- *Multiomic analysis figures for CD4 T cells in Muon*
- **Figure6_S12/Figure_6F_CD4na_DAPs.ipynb**
- *ATAC: Differentially Accessable Peaks in CD4 Naive T cells*

#### Figure7_S13_Treatment_effect
- **Figure7_S13_Treatment_effect/Figure_7BD_ABT_treatment.ipynb**
- *Abatacept responders and Non-responders analysis plots*
- **Figure7_S13_Treatment_effect/Figure_7CE_TNFi_Treatment.ipynb**
- *Abatacept TNFi treatment plots*
- **Figure7_S13_Treatment_effect/Figure_7H.ipynb**
- *Gene signature validation for Abatacept treatment*
- **Figure7_S13_Treatment_effect/Figure_S13_AB.ipynb**
- *DEGs in Abatacept Responders and Non-responders*

#### FigureS15
- **FigureS15/FigS15_A_B_C.ipynb**
- *Celltype score plot*
- **FigureS15/FigS15_D.ipynb**
- *Batch 182 effect in major celltypes*

## Usage

Clone the repository and run the notebooks in order. See `environment.yml` for dependencies.

For questions, contact the repository maintainer.
