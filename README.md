# scRNA-Seq Integration and Batch Effect Correction using Seurat

This repository contains an R script to integrate **single-cell RNA sequencing (scRNA-seq)** datasets and correct for **batch effects** using the [Seurat](https://satijalab.org/seurat/) package.  
The workflow is applied to hepatoblastoma datasets to identify distinct tumor cell populations and genetic mechanisms.

---

# Dataset Information

**Series:** [GSE180665](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180665)  
**Title:** Modeling Hepatoblastoma through Single Cell Sequencing (scRNA-seq)  
**Organism:** *Homo sapiens*  
**Type:** Expression profiling by high throughput sequencing  

---

## Summary
Hepatoblastoma (HB) is the most common childhood liver cancer. Using **patient-derived xenografts (PDX)**, this study explored tumor cell populations and pathways such as **WNT, AP-1, Hedgehog, Notch, and MAPK**, with genes like **GPC3, DLK1, and HMGA2** identified as key drivers. Distinct cell clusters linked to initiation, angiogenesis, maintenance, and progression were discovered, providing insights for novel therapeutic strategies.  

**Design:** Patient tumor implanted into **NSG mice**, comparing PDX tumor, background liver, and primary HB tumor with **scRNA-seq**.  
**Comparisons:** PDX tumor vs background liver vs primary HB tumor. 
**Dataset size:** ~26,886 cells across 3 samples after QC.  

---

## ⚙️ Workflow Overview

1. **Setup and Data Import**
   - Reads multiple 10X Genomics matrices.
   - Creates Seurat objects for each dataset.
   - **Merges all samples first**, allowing one-time QC across all datasets (efficient, avoids repetitive filtering).

2. **Quality Control (QC)**
   - Calculates **mitochondrial gene percentage (`percent.mt`)**.
   - Filters low-quality cells based on:
     - `nCount_RNA > 800`
     - `nFeature_RNA > 500`
     - `percent.mt < 10%`
   - Generates QC plots:
     - Distribution of RNA features, counts, and mitochondrial content.
     - Feature scatter plots.

3. **Pre-Integration Analysis**
   - Normalization, scaling, PCA, clustering, and UMAP.
   - Visualizations based on **Patient ID** and **Sample Type**.
   - Identifies batch effects (patient-driven clustering).

4. **Batch Effect Correction (Integration)**
   - Splits data by patient.
   - Normalizes and selects integration features.
   - **Seurat v5 fix (`JoinLayers`) applied** to resolve multi-layer assay errors.  
     This is important for others using Seurat v5, as `GetAssayData()` fails without it.
   - Finds integration anchors using **Canonical Correlation Analysis (CCA)**.
   - Integrates datasets into a unified space.

5. **Post-Integration Analysis**
   - Scaling, PCA, UMAP.
   - Visualization of clusters by **Patient ID** and **Sample Type**.
   - Side-by-side comparison of pre- and post-integration embeddings.

---

## 🖼️ Key Visualizations & Interpretations

---

### 1. QC Plots
<img width="652" height="543" alt="QC Plot" src="https://github.com/user-attachments/assets/963ffe5f-540d-4b98-9729-9b2de22a0939" />

The histogram shows a **bimodal distribution**, with:  
- A large peak around **500–1000 features** (likely low-quality cells or debris).  
- A second, broader peak between **2000–4000 features**.  

This indicates the necessity of filtering to remove cells in the lower peak.

---

### 2. Feature Scatter Plots
<img width="652" height="543" alt="Feature Scatter" src="https://github.com/user-attachments/assets/3a329245-8eec-411b-b7f0-0c3b47705ef0" />

- **Left Plot: nCount_RNA vs. nFeature_RNA**  
  A high positive correlation (**0.9**) is expected, as more unique transcripts (counts) should correlate with more unique genes (features).  
  The strong linear relationship confirms this.  

- **Right Plot: nCount_RNA vs. percent.mt**  
  A low/negative correlation is desired. The observed **-0.2 correlation** and the cluster of cells with low %mt indicates that cell complexity is not strongly driven by mitochondrial contamination.

---

### 3. Pre-Integration UMAP Interpretation  

#### 3.1 Elbow Plot  
<img width="651" height="544" alt="Elbow Plot" src="https://github.com/user-attachments/assets/27a119fe-d244-4024-86cf-3ef64d518cfb" />

Shows the standard deviation of principal components (PCs).  
The "elbow," where the drop-off in standard deviation sharply decreases, is typically around **PC 10–12**.  
This guided the choice of using the **top 20 PCs** for clustering and UMAP (a conservative, inclusive choice).

#### 3.2 UMAP (Pre-Integration)  
<img width="651" height="544" alt="UMAP Pre-Integration" src="https://github.com/user-attachments/assets/7eab33c9-7c37-4df8-a439-ab88714f008a" />

This plot shows the **38 clusters (0–37)** identified before integration.

#### 3.3 UMAP by Patient and Cell Type  
<img width="651" height="544" alt="UMAP Patient and Type" src="https://github.com/user-attachments/assets/1718171e-161b-4dab-8047-eb8e4449cd30" />

- **Patient (Left):** Cells predominantly cluster by patient (batch).  
  For example, cells from **HB30 (Green)** form a large, distinct group separate from **HB53 (Blue)**, while **HB17 (Red)** forms its own lower cluster.  
  This is the definition of a **batch effect**, where technical variation (patient/sample source) dominates true biological signal.  

- **Type (Right):** Although heavily influenced by batch effects, distinct biological populations (Tumor, PDX, Background) are discernible.  
  The **goal of integration** is to bring the patient-specific clusters closer together to enhance biological separation.

---

## 🔍 PCA Loadings: Before vs After Integration

Principal Component Analysis (PCA) was performed **before and after integration** to evaluate how batch effects influenced gene variation across patients.

### 📊 PCA Before Integration
- **PC1–PC3** are dominated by **liver metabolism and patient-specific genes**:  
  - High loadings of *CYP family* (CYP2B6, CYP3A4, CYP2C9, CYP3A5), *ABCB11*, *CP*, *CRP*, *APOB*, etc.  
  - Patient-driven signals (e.g., *GPC3, PEG10, AFP*) dominate negative axes.  
- **Interpretation:** PCA mainly captured **technical and sample-of-origin variation**, not just biological differences.  
  Clusters in UMAP separated strongly by patient → evidence of batch effect.

---

### 📊 PCA After Integration
- **PC1–PC2** are now dominated by **immune-related and developmental genes** rather than metabolic/patient-specific markers:  
  - Positive axes enriched for transcription factors and signaling (*NKD1, PDE4D, HMGA2, PTPRG, TCF4, MEF2C*).  
  - Negative axes enriched for housekeeping/liver function genes (*ALB, CPS1, FGB, FGG, APOC3*).  
- **PC3–PC5** capture tumor markers and lineage-specific variation:  
  - Tumor-associated (*AFP, PEG10, LGR5, CASC9*) now appear as biological drivers.  
  - Immune/inflammatory genes (*DOCK2, ZEB2, AOAH, CD86, MS4A6A*) indicate microenvironment contributions.  
- **Interpretation:** After integration, PCA reflects **biological variation (tumor vs microenvironment, differentiation states)** rather than patient-specific batch effects.

---

### 🧾 PCA Summary
- **Before integration:** PCA was biased by patient identity and batch-specific expression (metabolic genes, liver-specific background).  
- **After integration:** PCA highlights **true biology** — tumor markers (*AFP, PEG10, LGR5*), immune response genes (*ZEB2, CD86, DOCK2*), and signaling pathways, demonstrating effective batch correction.

---

### 4. Post-Integration UMAP Interpretation  

#### 4.1 UMAP after Integration by Patient and Cell Type  
<img width="701" height="504" alt="UMAP Post-Integration Patient and Type" src="https://github.com/user-attachments/assets/4a4f1605-600d-48b5-b1ef-1cf7a867d050" />

- **Patient (Left):** Batch effects are substantially corrected.  
  Cells from **HB17 (Red)**, **HB30 (Green)**, and **HB53 (Blue)** are now well-mixed within the major clusters.  
  This indicates that the corrected expression values now allow cells of the same biological type, regardless of the patient of origin, to group together.  

- **Type (Right):** The biological structure is preserved and enhanced.  
  The main clusters clearly separate into **Tumor (Maroon)**, **PDX (Purple)**, and **Background (Orange)** populations, confirming that the integration process successfully removed technical noise while retaining biological variance.

#### 4.2 Pre-Integration and Post-Integration by Patient and Cell Type (Comparison)  
<img width="701" height="504" alt="grid_after_intg" src="https://github.com/user-attachments/assets/8e61b02d-2268-4e6d-a432-bb89a68f5f79" />

This figure provides a direct comparison of the UMAP visualization before (Top Row) and after (Bottom Row) running the data integration steps, serving as the final justification for the entire Seurat integration pipeline.

#### 4.3 Final UMAP (Post-Integration)  
<img width="701" height="504" alt="Final UMAP" src="https://github.com/user-attachments/assets/43ca5be1-4230-483f-8d58-65f3dfb9fe97" />

This final UMAP shows the **refined clustering structure after integration**, resulting in **38 distinct clusters (0–37)**.  
The contiguous nature of these clusters confirms the high quality of the integrated data, which is now ready for **accurate cell type annotation** and subsequent biological analysis.

---

## 📦 Dependencies

- **R (≥ 4.1.0)**
- [Seurat](https://satijalab.org/seurat/) (v5.x recommended)
- `ggplot2`
- `tidyverse`
- `gridExtra`

Install missing packages with:

```R
install.packages(c("ggplot2", "tidyverse", "gridExtra"))
if (!requireNamespace("Seurat", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat")
}
