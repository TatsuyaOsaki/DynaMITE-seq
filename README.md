## R script for M. Campisi et al. Cancer Cell, 44, 1-21 (2026) 
Repository for the code that was used to process and analyze data collected for "Vascular STING activation facilitates NK cell anti-tumor immunity in small cell lung cancer", M. Campisi et al. Cancer Cell, 44, 1-21 (2026) 

------------------------------------------------------------------------

# scRNAseq analysis: DynaMITE-seq Single-Cell and Spatial Integration Pipeline

This repository contains a structured analysis workflow for processing
and integrating **DynaMITE-seq single-cell RNA-seq data** with spatial
transcriptomics and downstream pathway enrichment analysis.

# Workflow Overview

The pipeline is organized into three sequential Jupyter notebooks:

## 1. scRNAseq_DynaMITEseq.ipynb

### Purpose:
Preprocessing and analysis of DynaMITE-seq single-cell RNA-seq data.

### Main steps 
- Quality control (filtering low-quality cells and genes)
- Normalization and scaling
- Dimensionality reduction (PCA, UMAP)
- Clustering
- Cell type annotation
- Differential gene expression analysis

This notebook generates the processed dataset used for downstream
integration.

------------------------------------------------------------------------

## 2. Integrate to Visim.ipynb

### Purpose:
Integration of processed single-cell data with spatial transcriptomics
(Visium).

### Main steps 
- Data harmonization between scRNA-seq and spatial datasets
- Label transfer or anchor-based integration
- Mapping cell types to spatial coordinates
- Visualization of spatial gene expression patterns

This step links cellular transcriptomic resolution with spatial context.


## 3.ssGSEA_DynaMITEseq.ipynb

### Purpose:
Single-sample Gene Set Enrichment Analysis (ssGSEA).

### Main steps 
- Selection of curated gene sets (e.g., pathway databases)
- Calculation of ssGSEA enrichment scores
- Comparison of pathway activity across clusters or conditions
- Visualization of enrichment patterns

This notebook provides pathway-level functional interpretation.

### Pipeline Structure

Raw scRNA-seq\
→ QC & Clustering\
→ Spatial Integration\
→ Pathway Enrichment

------------------------------------------------------------------------
## Cell chat analysis

## 4. Cellchat_compute_interaction.ipynb

### Purpose
Inference of cell–cell communication networks using CellChat.

### Main Steps
- Construction of CellChat object from processed single-cell or spatial data  
- Identification of ligand–receptor interaction pairs  
- Calculation of communication probability scores  
- Aggregation of signaling pathways  
- Network visualization (interaction strength, outgoing/incoming signaling roles)  

### Output
- Cell–cell interaction matrices  
- Pathway-specific signaling networks  
- Network centrality metrics (sender/receiver/modulator roles)  

This notebook quantitatively models intercellular communication based on curated ligand–receptor databases.

---

## 5. Cellchat_differential_analysis.ipynb

### Purpose
Comparative analysis of cell–cell communication between biological conditions.

### Main Steps
- Merging CellChat objects across conditions (e.g., control vs disease)  
- Differential interaction strength analysis  
- Identification of altered signaling pathways  
- Comparative network topology assessment  
- Visualization of gained/lost signaling interactions  

### Output
- Differential ligand–receptor interactions  
- Condition-specific pathway activity shifts  
- Network-level rewiring metrics  

This notebook enables mechanistic interpretation of condition-dependent communication changes.

---

## Pipeline Structure

Raw scRNA-seq  
→ QC & Clustering (1)  
→ Spatial Integration (2)  
→ Pathway Enrichment (3)  
→ Cell–Cell Interaction Inference (4)  
→ Differential Communication Analysis (5)  

# Visium analysis

TBA
```
Marco Campisi1, Tatsuya Osaki2†, Ian Dryg3,4†, Carla Stornante1,5, Jacqueline Wolff3, Jason Weirather3,4, Nicholas Weaver3, Mubin Tarannum1, Ian Gillanders1, Alan Bers1, Eden Bobilev1, Maia Shea Lineberry1, Pieter Schol1#, Minyue Chen1#, Keiichi Ota1, Sung Rye Park6, Caroline G. Fahey1,7, Tran C. Thai1, Yixiang Li1, Zhaorong Li1, Ze-Hua Li1, Sarah Shelton8, Satomi Hirose9, Clare Padrick1,10, Elena Ivanova1,10, Minh Ha1,10, Shaobo Yang1, Harrison Olszewski3, Sophie Kivlehan1,10, Patrick Lizotte1,10, Kelsey Hagen1,10, Iliana Gjeci1,10, William Haller1,10, Lauren M. Zasadil1,10, Talal El Zarif1, Jacob E. Berchuck11, Erik H. Knelson1#, Elliott J. Brea1,12, Igor Odintsov13, Tianju Liu14, Valeria Chiono5, Bishma Tuladhar1,10, Kenneth Ngo1,10, Emily Rosentrater15, Jeffrey Raizer15, Neil Lineberry15, Mriganka Sur2, Kathleen Pfaff3, Matthew Freedman1, F. Stephen Hodi1,3, Michael Y. Tolstorukov6, Matthew G. Oser1, Michal Sheffer1, Constantine S. Mitsiades1,16,17,18, Michael J. Barnes19, Vicky A. Appleman15, Prafulla C. Gokhale1,10, Roger D. Kamm8,9, Cloud P. Paweletz1,10, Rizwan Romee1,20, Scott Rodig3,13, David A. Barbie1,10,12,16,20*, Navin R. Mahadevan1,13,14*

*Integrative profiling of primary and microphysiological tumor environments unveils a STING-mediated vascular checkpoint modulating NK cell immunity*
1.Department of Medical Oncology, Dana-Farber Cancer Institute, Boston, MA, 02215, USA
2.Picower Institute for Learning and Memory, Massachusetts Institute of Technology, Cambridge, MA, 02139, USA
3.Center for Immuno-Oncology, Dana-Farber Cancer Institute, Boston, MA, 02215, USA
4.Department of Data Science, Dana-Farber Cancer Institute, Boston, MA, 02215, USA
5.Department of Mechanical and Aerospace Engineering, Politecnico di Torino, Turin, 10129, Italy
6.Department of Informatics and Analytics, Dana-Farber Cancer Institute, Boston, MA, 02215, USA
7.Harvard University, Cambridge, MA, 02138, USA
8.Department of Biological Engineering, Massachusetts Institute of Technology, Cambridge, MA, 02139, USA
9.Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA, 02139, USA
1.0Belfer Center for Applied Cancer Science, Dana-Farber Cancer Institute, Boston, MA, 02215, USA
11.Department of Hematology and Medical Oncology, Emory University School of Medicine, Atlanta, GA, 30307, USA
12.Lowe Center for Thoracic Oncology, Dana-Farber Cancer Institute, Boston, MA, 02215, USA
13.Department of Pathology, Brigham and Women’s Hospital, Boston, MA, 02115, USA
14.Department of Pathology, University of Michigan, Ann Arbor, MI, 48109, USA
15.Takeda Development Center Americas, Inc. (TDCA), Lexington, MA, 02421, USA
16.Harvard Medical School, Boston, MA, 02215, USA
17.Broad Institute of MIT and Harvard, Cambridge, MA, 02142, USA
18.Ludwig Center at Harvard, Boston, MA, 02215, USA
19.Cancer Immunology and Cell Therapy Research, Bristol Myers Squibb, Seattle, WA 98109
20.Parker Institute for Cancer Immunotherapy

```
