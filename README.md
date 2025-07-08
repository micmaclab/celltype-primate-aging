# Cellular basis for cortical network aging in primates

This repository hosts data, code, and links to external resources used to generate the results of this work.

## Publication Link

[Preprint Link](https://doi.org/)

## Alignment of Data & Analyses Across Scales 

<p align="center">
  <img src="figs/overview.png" width="850">
</p>


## Project Code

- **Cross-species alignment of functional domains**  
  Usage: Map human functional networks onto macaque cortex using joint embedding  
  Code: [folder name]

- **Cortical similarity networks using Morphometric Inverse Divergence (MIND)**  
  Usage: Compute region-wise similarity based on multivariate MRI features  
  Code: [folder name]

- **Computation of network-based properties**  
  Usage: Derive total similarity strength across cortical regions and networks  
  Code: [folder name]

- **Statistical modeling of age effects**  
  Usage: Model age-related change in regional similarity strength using `lme4`  
  Code: [folder name]

- **3D single-cell transcriptomic atlas of macaque monkey**  
  Usage: Analyze region- and layer-specific cell type distributions  
  Code: [folder name]

- **Univariate associations with cell type specific abundances**  
  Usage: Correlate regional MRI features with cell-type abundance using null models  
  Code: [folder name]

- **Canonical Correlation Analysis**  
  Usage: Identify multivariate relationships between cell profiles and MRI-based MIND.   
  Code: [folder name]

- **Human specific cortical cell type enrichment dataset**  
  Usage: Integrate imputed human cell type distributions from Zhang et al (see *External Repositories*) 
  Code: [folder name]

- **Comparison of macaque and human cell types by their aligned functional networks**  
  Usage: Test for cell type enrichment within functional networks in both species  
  Code: [folder name]

- **Macaque to human comparison of cell type composition**  
  Usage: Compare cross-species cell type composition using PCA and Mantel test  
  Code: [folder name]



## External Repositories

- **MIND network construction**  
  - Code: https://github.com/isebenius/MIND  
  - Reference: [Sebenius, I. et al., *Nature Neuroscience* 26, 1461–1471 (2023).](https://doi.org/10.1038/s41593-023-01376-7)

- **Functional network cross-species alignment analyses**  
  - Code: https://github.com/TingsterX/alignment_macaque-human  
  - Reference: [Xu, T. et al., *NeuroImage* 223, 117346 (2020).](https://doi.org/10.1016/j.neuroimage.2020.117346)

- **3D single-cell transcriptomic atlas of macaque monkey**  
  - [Online Visualization Tool & Data](https://macaque.digital-brain.cn/spatial-omics)  
  - Reference: [Chen, A. et al., *Cell* 186, 3743.e24 (2023).](https://doi.org/10.1016/j.cell.2023.06.009)

- **Canonical correlation analysis (CCA)**  
  - Code: https://github.com/andersonwinkler/PermCCA  
  - Reference: [Winkler, A. M., et al., *NeuroImage* 220, 117065 (2020).](https://doi.org/10.1016/j.neuroimage.2020.117065)

- **Imputed human cell type distributions**  
  - Code: https://github.com/XihanZhang/human-cellular-func-con/tree/main
  - Reference: [Zhang et al., *Nature Neuroscience* 28, 150–160 (2025).](https://doi.org/10.1038/s41593-024-01812-2)

- **Surface Based Visualizations**  
  - Code: https://github.com/danjgale/surfplot  

## External Tools & Software Libraries

- **CIVET-macaque**  
  - Usage: Cortical surface extraction  
  - Reference: [Lepage, C. et al., *NeuroImage* 227, 117622 (2021).](https://doi.org/10.1016/j.neuroimage.2020.117622)

- **Connectome Workbench**  
  - Usage: Surface-based rendering and visualization  
  - Reference: [Marcus, D. et al., *Frontiers in Neuroinformatics* 5, (2011).](https://doi.org/10.3389/fninf.2011.00004)

- **netneurotools**  
  - Usage: Spatial permutation tests  
  - Language: Python  
  - Code: https://netneurolab.github.io/netneurotools/

- **nibabel**  
  - Usage: Neuroimaging data handling (NIfTI, GIFTI, etc.)  
  - Language: Python  
  - Code: https://nipy.org/nibabel/

- **neuromaps**  
  - Usage: Coordinate system registration and map comparison  
  - Language: Python  
  - Code: https://netneurolab.github.io/neuromaps/

- **lme4 (R package)**  
  - Usage: Linear mixed-effects modeling  
  - Language: R  
  - Code: https://cran.r-project.org/web/packages/lme4/index.html  
  - Reference: [Bates et al., *Journal of Statistical Software* 67(1), 2015.](https://doi.org/10.18637/jss.v067.i01)














