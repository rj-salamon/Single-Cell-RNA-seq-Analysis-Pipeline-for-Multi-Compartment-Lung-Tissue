# scRNA-seq Lung Analysis: Saline vs 6OHDA

## Repository overview

This repository contains scripts for single-cell RNA-seq (scRNA-seq) analysis of saline vs 6OHDA mouse lung datasets. The workflow includes preprocessing, integration, compartment-level annotation, mesenchymal and immune subclustering, mono/macrophage trajectory inference, and CellChat ligand–receptor communication analysis.

The analysis is implemented primarily using **Seurat**, **Slingshot**, **tradeSeq**, and **CellChat** in R.

---

## Data processing summary

Raw FASTQ files were processed using **10x Genomics Cell Ranger** against the mouse reference genome (**mm10**) to generate feature-barcode matrices. Downstream analysis was performed in R using **Seurat** for quality control, normalization, integration, clustering, and visualization.

Single-cell libraries were generated using the Chromium Single Cell 3’ v3 kit according to the manufacturer’s instructions. Sequencing was performed on an Illumina NovaSeq X Plus platform at the UC San Diego IGM Genomics Center.

---

## Directory structure

config/
data/
results/
figures/
scripts/
notebooks/
README.md


- **config/**  
  Configuration files (if applicable).

- **data/**  
  Optional location for feature-barcode matrices or intermediate data files not tracked by git.

- **results/**  
  Stores Seurat objects, differential expression tables, trajectory outputs, and CellChat results.

- **figures/**  
  Output figures (UMAPs, dot plots, volcano plots, heatmaps).

- **scripts/**  
  R scripts for each analysis step.

- **notebooks/**  
  Optional exploratory or visualization notebooks.

## Required inputs
The scripts assume the following input object exists in `results/`:

results/seurat_obj_mesenchymal_immune_merged_annotated.rds


This Seurat object contains four major compartments:
- Mesenchymal
- Epithelial
- Endothelial
- Immune

### Required metadata columns
- `compartment`
- `orig.ident` (e.g. `saline`, `6OHDA`)
- `mes_identity` (mesenchymal subtypes)
- `imm_identity` (immune major identities)
- `subcluster_identity` (mono/mac refined identities; `NA` for non-mono/mac immune cells)

If `subcluster_identity` is missing, the CellChat analysis will fall back to `imm_identity` for immune cells.

## Outputs
Each script writes outputs to `results/` and figures to `figures/` by default.

### Example outputs

**Reclustering / annotation**
- `results/seurat_obj_immune_subset.rds`
- `results/seurat_obj_mesenchymal_immune_merged_annotated.rds`

**Trajectory**
- `results/mono_slingshot_sce.rds`
- `results/tradeSeq_assoc_results.csv`

**CellChat**
- `results/cellchat_saline.rds`
- `results/cellchat_6OHDA.rds`
- `results/cellchat_merged.rds`
- `results/cellchat_top_pathways.RData`

## How to run
Run the analysis in order using the `.Rmd` files in `notebooks/` (or render them to HTML). Each step writes outputs to `results/` and figures to `figures/` by default.

### 1. Preprocess
`notebooks/01_scRNAseq_seurat_preprocess_rjs.Rmd`  
   - Loads raw feature-barcode matrices (e.g., 10x output)  
   - Performs initial QC filtering, normalization, and preprocessing  
   - Saves a processed Seurat object for downstream integration

### 2. Integration and clustering
`notebooks/02_integration_clustering_rjs.Rmd`  
   - Integrates samples/conditions (saline vs 6OHDA)  
   - Runs dimensional reduction (PCA/UMAP) and clustering  
   - Produces an integrated object used for compartment annotation/subsetting

### 3. Compartment level clustering
`notebooks/03_cluster_subset_4_compartment_rjs.Rmd`  
   - Defines four major compartments (Mesenchymal, Epithelial, Endothelial, Immune)  
   - Generates compartment-level UMAPs and summary plots  
   - Saves the 4-compartment annotated Seurat object used for lineage-specific analyses

### 4. Mesenchymal clustering
`notebooks/04_mesenchymal_subcluster_analysis_rjs.Rmd`  
   - Subsets mesenchymal cells  
   - Reclusters and annotates mesenchymal subtypes (`mes_identity`)  
   - Exports mesenchymal cluster markers and/or differential expression outputs

### 5. Immune clustering
`notebooks/05_immune_subcluster_analysis_rjs.Rmd`  
   - Subsets immune cells  
   - Reclusters and annotates immune identities (`imm_identity`)  
   - Exports immune cluster markers and/or differential expression outputs

### 6. Monocyte and Macrophage subclustering
`notebooks/06_mono-mac_immune_subcluster_rjs.Rmd`  
   - Subsets mono/mac populations from the immune compartment  
   - Performs higher-resolution mono/mac subclustering  
   - Writes refined labels back into the merged object as `subcluster_identity`

### 7. Trajectory analysis in immmune 
`notebooks/07_trajectory_analysis_immune_rjs.Rmd`  
   - Runs Slingshot pseudotime analysis on mono/mac subclusters  
   - Fits GAMs with tradeSeq to identify pseudotime-associated genes  
   - Saves trajectory objects and result tables/figures

### 8. Cell chat mesenchymal-immune
`notebooks/08_cell_chat_immune_mesenchymal_rjs.Rmd`  
   - Builds CellChat objects per condition (saline vs 6OHDA)  
   - Uses a unified `celltype_fine` label:
     - Mesenchymal: `mes_identity`
     - Immune: `subcluster_identity` when available, otherwise `imm_identity`
   - Identifies top signaling pathways and differential communication between conditions

## Key parameters

### Seurat clustering
- SCT normalization (`SCTransform`, `glmGamPoi` where specified)

- PCA dimensions:
  - Mesenchymal reclustering: typically 1:30
  - Immune reclustering: typically 1:30
  - Mono/mac subclustering: typically 1:20

- Clustering resolution:
  - Mesenchymal: ~0.4
  - Immune: ~0.4
  - Mono/mac: ~0.15–0.2

### Differential expression filtering
- `min.pct = 0.25`
- `FDR < 0.05`
- `|avg_log2FC| > 0.5` (biological and visualization threshold)

### CellChat
- `min.cells = 10` in `filterCommunication`

## Reproducibility
R package versions can be captured using either:
- `sessionInfo()` (saved to `results/sessionInfo.txt`), or  
- `renv` (recommended), generating an `renv.lock` file

## Notes on cell type labels used for CellChat
CellChat grouping is defined using a single metadata column:

- `celltype_fine`


- Mesenchymal cells: `mes_identity`
- Immune cells: `subcluster_identity` when available (mono/mac refined labels),
  otherwise `imm_identity`

This strategy preserves the full immune atlas while enabling higher-resolution
communication analysis for mono/macrophage subtypes.

---

## Reproducibility

R package versions can be captured using either:

- `sessionInfo()` (saved to `results/sessionInfo.txt`), or  
- `renv` (recommended), generating an `renv.lock` file

