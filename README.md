# CAZyme Fold Architecture Is Conserved Between Disparate Environments Despite Extreme Sequence Divergence

**Accepted at *mSystems* (2026)**
DOI: [10.1128/msystems.00485-26](https://doi.org/10.1128/msystems.00485-26) (in press)

Oliyad Jeilu, Addis Simachew, Erica M. Hartmann, Erik Alexandersson, Eva Johansson

---

## Overview

This study compares the carbohydrate-degradation potential of two contrasting ecosystems using shotgun metagenomics: the alkaline-saline soda lakes of the East African Rift Valley (Lakes Abijata, Chitu, and Shala) and the anaerobic ruminant gut (goat, cattle, and sheep). From 34 metagenomes, we recovered 371 quality-filtered MAGs, annotated 26,541 CAZyme genes, and predicted three-dimensional structures for 12 representative glycoside hydrolases using AlphaFold3. We show that despite extreme sequence divergence between environments, the catalytic fold architecture of key CAZymes is remarkably conserved.

---

## Key Findings

- Soda lake MAGs harbor greater phylogenetic novelty (84% novel species vs. 52% in rumen) and lower RED values, indicating deeper evolutionary divergence from classified reference genomes.
- Rumen MAGs are enriched in fibrolytic GH families (GH2, GH9, GH10, GH28), while soda lake MAGs are enriched in cell wall remodeling enzymes (GH23, GH73, GH103) and stress-response pathways (ectoine biosynthesis, Sox sulfur oxidation, nitrate/nitrous oxide reduction).
- All 12 AlphaFold3-predicted structures adopted canonical GH family folds with high confidence (pTM 0.75--0.97), and 23 of 26 catalytic residues were conserved across both environments.
- Soda lakes and rumen metagenomes represent complementary reservoirs of industrially relevant CAZyme diversity.

---

## Repository Structure

```
.
├── README.md
├── soda_rumen_metagenomics_analysis.Rmd   # Complete reproducible analysis (published version)
├── Bash_scripts_upstream/                  # HPC pipeline scripts
│   ├── QC.sh                               # Read quality control (FastQC, Trimmomatic, KneadData)
│   ├── Assembly_Megahit.sh                 # Metagenomic assembly (MEGAHIT)
│   ├── Quast.sh                            # Assembly quality assessment (QUAST)
│   ├── binning_MEtabat.sh                  # MAG binning (MetaBAT2)
│   ├── MAGs taxonomy_GTDBk.sh              # Taxonomic classification (GTDB-Tk v2.4.0)
│   ├── Phylogenetic tree_Phylophlan.sh     # Phylogenomic tree (PhyloPhlAn v3.0)
│   ├── Functional annotation_Prokka.sh     # Gene prediction (Prokka v1.14.5)
│   ├── CAZy annoattion.sh                  # CAZyme annotation (dbCAN2)
│   ├── taxonmoy_metaphlan.sh               # Taxonomic profiling (MetaPhlAn4)
│   ├── function_humman3.sh                 # Functional profiling (HUMAnN3)
│   └── Function prediction based on structure_ECOFOLD.sh  # EcoFoldDB annotation
│
└── R scripts Downstream analysis/          # R scripts (earlier versions)
    ├── taxonomy_abundance_alpha_beta_analysis.R
    ├── Mags stat taxonomy.R
    ├── Cazy_ecofold_analysisi.R
    ├── humman3 metacyc.R
    ├── phylogentictree.R
    └── R_scripts/
```

---

## Main Analysis Pipeline

The complete reproducible analysis is in **`soda_rumen_metagenomics_analysis.Rmd`**, which generates all 7 main figures and Table 1 from the manuscript:

| Figure | Description |
|--------|-------------|
| Fig 1 | Community composition: phylum bars, genus heatmap, unclassified fraction, alpha/beta diversity |
| Fig 2 | MAG quality, RED, ANI, taxonomic novelty, phylum/species composition |
| Fig 3 | Circular phylogenomic tree of 245 dereplicated MAGs |
| Fig 4 | Functional pathway prevalence heatmap (29 pathways, 9 categories) |
| Fig 5 | CAZyme class distribution, per-MAG diversity, sequence identity |
| Fig 6 | GH family repertoire heatmap (15 GH families across top genera) |
| Fig 7 | AlphaFold3-predicted structures for 12 CAZymes from 6 GH families |

---

## Data Availability

Raw sequencing data: NCBI SRA BioProject [PRJNA1273195](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1273195)

---

## Analysis Environment

Upstream bioinformatics was performed on Swedish national HPC clusters:

- **Dardel**, PDC Center for High Performance Computing (KTH Royal Institute of Technology)
- **UPPMAX**, Uppsala Multidisciplinary Center for Advanced Computational Science

Downstream statistical analyses and visualizations were conducted in R v4.3.1 using vegan, ComplexHeatmap, ggplot2, ggtree, and patchwork. Protein structures were predicted with AlphaFold3 and compared using Foldseek against PDB100.

---

## Citation

If you use this code or data, please cite:

> Jeilu O, Simachew A, Hartmann EM, Alexandersson E, Johansson E. CAZyme fold architecture is conserved between disparate environments despite extreme sequence divergence. *mSystems*. 2026 (in press). DOI: [10.1128/msystems.00485-26](https://doi.org/10.1128/msystems.00485-26)

---

## License

This project is shared for academic use. Please cite the publication if you reuse any part of this work.
