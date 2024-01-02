# diHMM (C++/Python version)

## Note: 
This is the updated version of the diHMM model (Marco et al. 2017, Nat Comm). The original package was written in MATLAB and can be accessed at gcyuan/diHMM. In this updated version, we have increased computational efficiency in two major ways: 1. Implementing the code in C++ along with a Python wrapper to increase computational speed. 2. To use an ensemble approach to aggregate information from multiple models each trained on a different sample. The details of these changes are described in (Kai et al. 2020 bioRxiv) 

## Overview
diHMM stands for Hierarchical Hidden Markov Model. diHMM is a computational method for finding chromatin states at multiple scales. The model takes as input a multidimensional set of histone modifications for several cell types and classifies the genome into a preselected number of nucleosome-level and domain-level hidden states. 

The diHMM model was originally developed by Eugenio Marco, with assistance from Wouter Meuleman, Jialiang Huang, Luca Pinello, Manolis Kellis and Guo-Cheng Yuan. The method was originally implemented in MATLAB. The code and sample data are available at (https://github.com/gcyuan/diHMM).

In this updated version, the computational efficiency is siginificantly improved by implementing in C++. Additional improvement is achieved by using an ensemble clustering approach. The development of this newer version was led by Stephanos Tsoucas and Yan Kai with assistance from Shengbao Suo, Xuan Cao, and Guo-Cheng Yuan.

**References:**
Marco E*, Meuleman W*, Huang J*, Glass K, Pinello L, Wang J, Kellis M†, Yuan GC†. Multi-scale chromatin state annotation using a hierarchical hidden Markov model. Nature Commun. 2017 Apr 7;8:15011. (https://www.nature.com/articles/ncomms15011)).

Kai Y, Tsoucas S, Suo S, Yuan GC. Multi-scale annotations of chromatin states in 127 human cell-types. bioRxiv. (https://www.biorxiv.org/content/10.1101/2020.12.22.424078v1).


## Annotation of 127 human cell types

We applied diHMM v1.0 to generate the multi-scale chromatin state annotations for the 127 human reference epigenomes in the Roadmap and ENCODE consortia. Detailed information of the information on the 127 epigenomes can be found at [Roadmap Epigenomics Consortium and Kundaje et.al](https://www.nature.com/articles/nature14248) or [this site](https://egg2.wustl.edu/roadmap/web_portal/meta.html).

## A first impression of diHMM states
We generated the chromatin state maps at the nucleosome (200bp resolution) and domain (4kb resolution) level. Here we show a snapshot of diHMM states across 127 epigenomes below.

![example](images/genome-browser_final_for_github.png)

## Accessing the multi-scale chromatin state maps in the 127 epigenomes
Our multi-scale annotations can be freely downloaded from the table below. After unzipping, those maps can be directly uploaded to genome browsers (e.g IGV) for visualization. Please note that our state annotations are based on hg19 reference genome.

To see full meta information about each reference epigenome, please visit [here](https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15).

| Epigenome   ID (EID) | Nucleosome                                                                                                  | Domain                                                                                                  |        GROUP     |        Standardized Epigenome name                               |        ANATOMY     |
|----------------------|-------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|------------------|------------------------------------------------------------------|--------------------|
| E017                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E017_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E017_domain.bed.gz) | IMR90            | IMR90   fetal lung fibroblasts Cell Line                         | LUNG               |
| E002                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E002_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E002_domain.bed.gz) | ESC              | ES-WA7   Cells                                                   | ESC                |
| E008                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E008_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E008_domain.bed.gz) | ESC              | H9   Cells                                                       | ESC                |
| E001                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E001_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E001_domain.bed.gz) | ESC              | ES-I3   Cells                                                    | ESC                |
| E015                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E015_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E015_domain.bed.gz) | ESC              | HUES6   Cells                                                    | ESC                |
| E014                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E014_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E014_domain.bed.gz) | ESC              | HUES48   Cells                                                   | ESC                |
| E016                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E016_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E016_domain.bed.gz) | ESC              | HUES64   Cells                                                   | ESC                |
| E003                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E003_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E003_domain.bed.gz) | ESC              | H1   Cells                                                       | ESC                |
| E024                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E024_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E024_domain.bed.gz) | ESC              | ES-UCSF4   Cells                                                 | ESC                |
| E020                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E020_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E020_domain.bed.gz) | iPSC             | iPS-20b   Cells                                                  | IPSC               |
| E019                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E019_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E019_domain.bed.gz) | iPSC             | iPS-18   Cells                                                   | IPSC               |
| E018                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E018_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E018_domain.bed.gz) | iPSC             | iPS-15b   Cells                                                  | IPSC               |
| E021                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E021_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E021_domain.bed.gz) | iPSC             | iPS   DF 6.9 Cells                                               | IPSC               |
| E022                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E022_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E022_domain.bed.gz) | iPSC             | iPS   DF 19.11 Cells                                             | IPSC               |
| E007                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E007_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E007_domain.bed.gz) | ES-deriv         | H1   Derived Neuronal Progenitor Cultured Cells                  | ESC_DERIVED        |
| E009                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E009_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E009_domain.bed.gz) | ES-deriv         | H9   Derived Neuronal Progenitor Cultured Cells                  | ESC_DERIVED        |
| E010                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E010_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E010_domain.bed.gz) | ES-deriv         | H9   Derived Neuron Cultured Cells                               | ESC_DERIVED        |
| E013                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E013_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E013_domain.bed.gz) | ES-deriv         | hESC   Derived CD56+ Mesoderm Cultured Cells                     | ESC_DERIVED        |
| E012                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E012_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E012_domain.bed.gz) | ES-deriv         | hESC   Derived CD56+ Ectoderm Cultured Cells                     | ESC_DERIVED        |
| E011                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E011_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E011_domain.bed.gz) | ES-deriv         | hESC   Derived CD184+ Endoderm Cultured Cells                    | ESC_DERIVED        |
| E004                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E004_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E004_domain.bed.gz) | ES-deriv         | H1   BMP4 Derived Mesendoderm Cultured Cells                     | ESC_DERIVED        |
| E005                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E005_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E005_domain.bed.gz) | ES-deriv         | H1   BMP4 Derived Trophoblast Cultured Cells                     | ESC_DERIVED        |
| E006                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E006_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E006_domain.bed.gz) | ES-deriv         | H1   Derived Mesenchymal Stem Cells                              | ESC_DERIVED        |
| E062                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E062_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E062_domain.bed.gz) | Blood   & T-cell | Primary   mononuclear cells from peripheral blood                | BLOOD              |
| E034                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E034_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E034_domain.bed.gz) | Blood   & T-cell | Primary   T cells from peripheral blood                          | BLOOD              |
| E045                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E045_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E045_domain.bed.gz) | Blood   & T-cell | Primary   T cells effector/memory enriched from peripheral blood | BLOOD              |
| E033                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E033_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E033_domain.bed.gz) | Blood   & T-cell | Primary   T cells from cord blood                                | BLOOD              |
| E044                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E044_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E044_domain.bed.gz) | Blood   & T-cell | Primary   T regulatory cells from peripheral blood               | BLOOD              |
| E043                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E043_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E043_domain.bed.gz) | Blood   & T-cell | Primary   T helper cells from peripheral blood                   | BLOOD              |
| E039                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E039_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E039_domain.bed.gz) | Blood   & T-cell | Primary   T helper naive cells from peripheral blood             | BLOOD              |
| E041                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E041_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E041_domain.bed.gz) | Blood   & T-cell | Primary   T helper cells PMA-I stimulated                        | BLOOD              |
| E042                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E042_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E042_domain.bed.gz) | Blood   & T-cell | Primary   T helper 17 cells PMA-I stimulated                     | BLOOD              |
| E040                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E040_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E040_domain.bed.gz) | Blood   & T-cell | Primary   T helper memory cells from peripheral blood 1          | BLOOD              |
| E037                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E037_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E037_domain.bed.gz) | Blood   & T-cell | Primary   T helper memory cells from peripheral blood 2          | BLOOD              |
| E048                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E048_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E048_domain.bed.gz) | Blood   & T-cell | Primary   T CD8+ memory cells from peripheral blood              | BLOOD              |
| E038                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E038_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E038_domain.bed.gz) | Blood   & T-cell | Primary   T helper naive cells from peripheral blood             | BLOOD              |
| E047                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E047_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E047_domain.bed.gz) | Blood   & T-cell | Primary   T CD8+ naive cells from peripheral blood               | BLOOD              |
| E029                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E029_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E029_domain.bed.gz) | HSC   & B-cell   | Primary   monocytes from peripheral blood                        | BLOOD              |
| E031                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E031_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E031_domain.bed.gz) | HSC   & B-cell   | Primary   B cells from cord blood                                | BLOOD              |
| E035                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E035_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E035_domain.bed.gz) | HSC   & B-cell   | Primary   hematopoietic stem cells                               | BLOOD              |
| E051                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E051_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E051_domain.bed.gz) | HSC   & B-cell   | Primary   hematopoietic stem cells G-CSF-mobilized Male          | BLOOD              |
| E050                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E050_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E050_domain.bed.gz) | HSC   & B-cell   | Primary   hematopoietic stem cells G-CSF-mobilized Female        | BLOOD              |
| E036                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E036_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E036_domain.bed.gz) | HSC   & B-cell   | Primary   hematopoietic stem cells short term culture            | BLOOD              |
| E032                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E032_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E032_domain.bed.gz) | HSC   & B-cell   | Primary   B cells from peripheral blood                          | BLOOD              |
| E046                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E046_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E046_domain.bed.gz) | HSC   & B-cell   | Primary   Natural Killer cells from peripheral blood             | BLOOD              |
| E030                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E030_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E030_domain.bed.gz) | HSC   & B-cell   | Primary   neutrophils from peripheral blood                      | BLOOD              |
| E026                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E026_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E026_domain.bed.gz) | Mesench          | Bone   Marrow Derived Cultured Mesenchymal Stem Cells            | STROMAL_CONNECTIVE |
| E049                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E049_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E049_domain.bed.gz) | Mesench          | Mesenchymal   Stem Cell Derived Chondrocyte Cultured Cells       | STROMAL_CONNECTIVE |
| E025                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E025_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E025_domain.bed.gz) | Mesench          | Adipose   Derived Mesenchymal Stem Cell Cultured Cells           | FAT                |
| E023                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E023_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E023_domain.bed.gz) | Mesench          | Mesenchymal   Stem Cell Derived Adipocyte Cultured Cells         | FAT                |
| E052                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E052_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E052_domain.bed.gz) | Myosat           | Muscle   Satellite Cultured Cells                                | MUSCLE             |
| E055                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E055_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E055_domain.bed.gz) | Epithelial       | Foreskin   Fibroblast Primary Cells skin01                       | SKIN               |
| E056                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E056_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E056_domain.bed.gz) | Epithelial       | Foreskin   Fibroblast Primary Cells skin02                       | SKIN               |
| E059                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E059_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E059_domain.bed.gz) | Epithelial       | Foreskin   Melanocyte Primary Cells skin01                       | SKIN               |
| E061                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E061_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E061_domain.bed.gz) | Epithelial       | Foreskin   Melanocyte Primary Cells skin03                       | SKIN               |
| E057                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E057_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E057_domain.bed.gz) | Epithelial       | Foreskin   Keratinocyte Primary Cells skin02                     | SKIN               |
| E058                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E058_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E058_domain.bed.gz) | Epithelial       | Foreskin   Keratinocyte Primary Cells skin03                     | SKIN               |
| E028                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E028_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E028_domain.bed.gz) | Epithelial       | Breast   variant Human Mammary Epithelial Cells (vHMEC)          | BREAST             |
| E027                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E027_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E027_domain.bed.gz) | Epithelial       | Breast   Myoepithelial Primary Cells                             | BREAST             |
| E054                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E054_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E054_domain.bed.gz) | Neurosph         | Ganglion   Eminence derived primary cultured neurospheres        | BRAIN              |
| E053                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E053_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E053_domain.bed.gz) | Neurosph         | Cortex   derived primary cultured neurospheres                   | BRAIN              |
| E112                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E112_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E112_domain.bed.gz) | Thymus           | Thymus                                                           | THYMUS             |
| E093                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E093_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E093_domain.bed.gz) | Thymus           | Fetal   Thymus                                                   | THYMUS             |
| E071                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E071_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E071_domain.bed.gz) | Brain            | Brain   Hippocampus Middle                                       | BRAIN              |
| E074                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E074_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E074_domain.bed.gz) | Brain            | Brain   Substantia Nigra                                         | BRAIN              |
| E068                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E068_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E068_domain.bed.gz) | Brain            | Brain   Anterior Caudate                                         | BRAIN              |
| E069                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E069_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E069_domain.bed.gz) | Brain            | Brain   Cingulate Gyrus                                          | BRAIN              |
| E072                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E072_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E072_domain.bed.gz) | Brain            | Brain   Inferior Temporal Lobe                                   | BRAIN              |
| E067                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E067_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E067_domain.bed.gz) | Brain            | Brain   Angular Gyrus                                            | BRAIN              |
| E073                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E073_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E073_domain.bed.gz) | Brain            | Brain_Dorsolateral_Prefrontal_Cortex                             | BRAIN              |
| E070                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E070_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E070_domain.bed.gz) | Brain            | Brain   Germinal Matrix                                          | BRAIN              |
| E082                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E082_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E082_domain.bed.gz) | Brain            | Fetal   Brain Female                                             | BRAIN              |
| E081                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E081_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E081_domain.bed.gz) | Brain            | Fetal   Brain Male                                               | BRAIN              |
| E063                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E063_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E063_domain.bed.gz) | Adipose          | Adipose   Nuclei                                                 | FAT                |
| E100                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E100_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E100_domain.bed.gz) | Muscle           | Psoas   Muscle                                                   | MUSCLE             |
| E108                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E108_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E108_domain.bed.gz) | Muscle           | Skeletal   Muscle Female                                         | MUSCLE             |
| E107                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E107_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E107_domain.bed.gz) | Muscle           | Skeletal   Muscle Male                                           | MUSCLE             |
| E089                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E089_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E089_domain.bed.gz) | Muscle           | Fetal   Muscle Trunk                                             | MUSCLE             |
| E090                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E090_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E090_domain.bed.gz) | Muscle           | Fetal   Muscle Leg                                               | MUSCLE_LEG         |
| E083                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E083_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E083_domain.bed.gz) | Heart            | Fetal   Heart                                                    | HEART              |
| E104                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E104_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E104_domain.bed.gz) | Heart            | Right   Atrium                                                   | HEART              |
| E095                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E095_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E095_domain.bed.gz) | Heart            | Left   Ventricle                                                 | HEART              |
| E105                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E105_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E105_domain.bed.gz) | Heart            | Right   Ventricle                                                | HEART              |
| E065                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E065_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E065_domain.bed.gz) | Heart            | Aorta                                                            | VASCULAR           |
| E078                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E078_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E078_domain.bed.gz) | Sm.   Muscle     | Duodenum   Smooth Muscle                                         | GI_DUODENUM        |
| E076                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E076_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E076_domain.bed.gz) | Sm.   Muscle     | Colon   Smooth Muscle                                            | GI_COLON           |
| E103                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E103_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E103_domain.bed.gz) | Sm.   Muscle     | Rectal   Smooth Muscle                                           | GI_RECTUM          |
| E111                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E111_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E111_domain.bed.gz) | Sm.   Muscle     | Stomach   Smooth Muscle                                          | GI_STOMACH         |
| E092                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E092_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E092_domain.bed.gz) | Digestive        | Fetal   Stomach                                                  | GI_STOMACH         |
| E085                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E085_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E085_domain.bed.gz) | Digestive        | Fetal   Intestine Small                                          | GI_INTESTINE       |
| E084                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E084_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E084_domain.bed.gz) | Digestive        | Fetal   Intestine Large                                          | GI_INTESTINE       |
| E109                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E109_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E109_domain.bed.gz) | Digestive        | Small   Intestine                                                | GI_INTESTINE       |
| E106                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E106_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E106_domain.bed.gz) | Digestive        | Sigmoid   Colon                                                  | GI_COLON           |
| E075                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E075_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E075_domain.bed.gz) | Digestive        | Colonic   Mucosa                                                 | GI_COLON           |
| E101                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E101_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E101_domain.bed.gz) | Digestive        | Rectal   Mucosa Donor 29                                         | GI_RECTUM          |
| E102                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E102_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E102_domain.bed.gz) | Digestive        | Rectal   Mucosa Donor 31                                         | GI_RECTUM          |
| E110                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E110_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E110_domain.bed.gz) | Digestive        | Stomach   Mucosa                                                 | GI_STOMACH         |
| E077                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E077_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E077_domain.bed.gz) | Digestive        | Duodenum   Mucosa                                                | GI_DUODENUM        |
| E079                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E079_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E079_domain.bed.gz) | Digestive        | Esophagus                                                        | GI_ESOPHAGUS       |
| E094                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E094_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E094_domain.bed.gz) | Digestive        | Gastric                                                          | GI_STOMACH         |
| E099                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E099_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E099_domain.bed.gz) | Other            | Placenta   Amnion                                                | PLACENTA           |
| E086                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E086_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E086_domain.bed.gz) | Other            | Fetal   Kidney                                                   | KIDNEY             |
| E088                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E088_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E088_domain.bed.gz) | Other            | Fetal   Lung                                                     | LUNG               |
| E097                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E097_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E097_domain.bed.gz) | Other            | Ovary                                                            | OVARY              |
| E087                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E087_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E087_domain.bed.gz) | Other            | Pancreatic   Islets                                              | PANCREAS           |
| E080                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E080_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E080_domain.bed.gz) | Other            | Fetal   Adrenal Gland                                            | ADRENAL            |
| E091                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E091_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E091_domain.bed.gz) | Other            | Placenta                                                         | PLACENTA           |
| E066                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E066_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E066_domain.bed.gz) | Other            | Liver                                                            | LIVER              |
| E098                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E098_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E098_domain.bed.gz) | Other            | Pancreas                                                         | PANCREAS           |
| E096                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E096_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E096_domain.bed.gz) | Other            | Lung                                                             | LUNG               |
| E113                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E113_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E113_domain.bed.gz) | Other            | Spleen                                                           | SPLEEN             |
| E114                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E114_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E114_domain.bed.gz) | ENCODE2012       | A549   EtOH 0.02pct Lung Carcinoma Cell Line                     | LUNG               |
| E115                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E115_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E115_domain.bed.gz) | ENCODE2012       | Dnd41   TCell Leukemia Cell Line                                 | BLOOD              |
| E116                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E116_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E116_domain.bed.gz) | ENCODE2012       | GM12878   Lymphoblastoid Cells                                   | BLOOD              |
| E117                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E117_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E117_domain.bed.gz) | ENCODE2012       | HeLa-S3   Cervical Carcinoma Cell Line                           | CERVIX             |
| E118                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E118_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E118_domain.bed.gz) | ENCODE2012       | HepG2   Hepatocellular Carcinoma Cell Line                       | LIVER              |
| E119                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E119_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E119_domain.bed.gz) | ENCODE2012       | HMEC   Mammary Epithelial Primary Cells                          | BREAST             |
| E120                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E120_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E120_domain.bed.gz) | ENCODE2012       | HSMM   Skeletal Muscle Myoblasts Cells                           | MUSCLE             |
| E121                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E121_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E121_domain.bed.gz) | ENCODE2012       | HSMM   cell derived Skeletal Muscle Myotubes Cells               | MUSCLE             |
| E122                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E122_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E122_domain.bed.gz) | ENCODE2012       | HUVEC   Umbilical Vein Endothelial Primary Cells                 | VASCULAR           |
| E123                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E123_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E123_domain.bed.gz) | ENCODE2012       | K562   Leukemia Cells                                            | BLOOD              |
| E124                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E124_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E124_domain.bed.gz) | ENCODE2012       | Monocytes-CD14+   RO01746 Primary Cells                          | BLOOD              |
| E125                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E125_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E125_domain.bed.gz) | ENCODE2012       | NH-A   Astrocytes Primary Cells                                  | BRAIN              |
| E126                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E126_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E126_domain.bed.gz) | ENCODE2012       | NHDF-Ad   Adult Dermal Fibroblast Primary Cells                  | SKIN               |
| E127                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E127_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E127_domain.bed.gz) | ENCODE2012       | NHEK-Epidermal   Keratinocyte Primary Cells                      | SKIN               |
| E128                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E128_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E128_domain.bed.gz) | ENCODE2012       | NHLF   Lung Fibroblast Primary Cells                             | LUNG               |
| E129                 | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E129_nucleosome.bed.gz) | [download](https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/E129_domain.bed.gz) | ENCODE2012       | Osteoblast   Primary Cells                                       | BONE               |

These maps can be freely downloaded from [here](https://www.dropbox.com/sh/85nxvu1hiwhwm9r/AAB0pQFvwD1KRqpwOOHf6A_Xa?dl=0).

## Installation
1. Create Conda Environment
   python version: 2.7
```
conda create -y -n dihmm python=2.7
conda activate dihmm
```

2. Downlaod dihmm-cpp
```
git clone https://github.com/gcyuan/diHMM-cpp.git
```

3. Install
Go into the build dir and run
```
cd diHMM-cpp/build
cmake ..
```
```
-- The C compiler identification is GNU 9.4.0
-- The CXX compiler identification is GNU 9.4.0
-- Detecting C compiler ABI info-- Detecting C compiler ABI info - done
-- Check for working C compiler: /sc/arion/projects/YuanLab/gcproj/xuan/anaconda3/envs/dihmm/bin/x86_64-conda-linux-gnu-cc - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /sc/arion/projects/YuanLab/gcproj/xuan/anaconda3/envs/dihmm/bin/x86_64-conda-linux-gnu-c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Found Armadillo: /sc/arion/projects/YuanLab/gcproj/xuan/anaconda3/envs/dihmm/include (found version "11.2.0") 
-- Found PythonLibs: /sc/arion/projects/YuanLab/gcproj/xuan/anaconda3/envs/dihmm/lib/libpython2.7.so (found suitable version "2.7.18", minimum required is "2.7") 
-- Found Boost: /sc/arion/projects/YuanLab/gcproj/xuan/anaconda3/envs/dihmm/lib/cmake/Boost-1.72.0/BoostConfig.cmake (found version "1.72.0") found components: python numpy filesystem 
-- Looking for sgemm_
-- Looking for sgemm_ - found
-- Found BLAS: /sc/arion/projects/YuanLab/gcproj/xuan/anaconda3/envs/dihmm/lib/libopenblas.so  
-- Configuring done
-- Generating done
-- Build files have been written to: /sc/arion/projects/YuanLab/gcproj/xuan/dihmm-cpp/build
```
```
make
```

```
Consolidate compiler generated dependencies of target dihmm
[ 14%] Building CXX object CMakeFiles/dihmm.dir/Model.cpp.o
[ 28%] Building CXX object CMakeFiles/dihmm.dir/Emissions.cpp.o
[ 42%] Building CXX object CMakeFiles/dihmm.dir/Forward_Backward.cpp.o
[ 57%] Linking CXX shared library libdihmm.so
[ 71%] Built target dihmm
Consolidate compiler generated dependencies of target dihmm_ext
[ 85%] Building CXX object CMakeFiles/dihmm_ext.dir/dihmm_ext.cpp.o
[100%] Linking CXX shared library dihmm_ext.so
[100%] Built target dihmm_ext
```

4. Install the dependency in your environment

bedtools, wigToBigWig, fetchChromSizes, bigWigToBedGraph required
   
```
conda install -c bioconda bedtools
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/fetchChromSizes
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigToBedGraph

```



5. Set path

```
export PYTHONPATH=${your_dihmm_dir}/diHMM-cpp/build
```

Then you can open a Python shell in the same dir and do
```
>>> import dihmm_ext
```

## Training a model
Training a diHMM model can be done by using the script *Train_diHMM.py*, after making necessary changes to input data path and other parameters.

Including applying a diHMM model for chromatin state annotation. Annotation can be done with the script *annotation.py* using in *Train_diHMM.py*.

```
python dihmm-cpp/Train_diHMM.py -h
usage: Train_diHMM.py [-h] -i INPUT_DIR --clusters CLUSTERS --chroms CHROMS
                         -o OUT_DIR [--n_bin_states N_BIN_STATES]
                         [--n_domain_states N_DOMAIN_STATES]
                         [--domain_size DOMAIN_SIZE] [--tolerance TOLERANCE]
                         [--max_iter MAX_ITER] [--bin_res BIN_RES]

Train diHMM runner.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR          The input binarized files dir. File name:
                        X1_chr1_binary.txt.
  --clusters CLUSTERS   Clusters/cell_types names used to train model.
                        Example: X1,X2 .
  --chroms CHROMS       chrs used to train model. Example: chr1,chr2 .
  -o OUT_DIR            Output dir.
  --n_bin_states N_BIN_STATES
                        Number of bin states. Default=2.
  --n_domain_states N_DOMAIN_STATES
                        Number of domain states. Default=4.
  --domain_size DOMAIN_SIZE
                        Number of domain size. Default=8.
  --tolerance TOLERANCE
                        Number of bin states. Default=1e-6.
  --max_iter MAX_ITER   Max iter number. Default=500.
  --bin_res BIN_RES     bin length used to generate binarized files.
                        Default=500.
```

## Tutorial
Here is the [tutorial](https://github.com/gcyuan/diHMM-cpp/blob/master/Tutorial/Example%20-%20Analyzing%20ENCODE%20H3K4me3%20data%20in%20hESC%20Cells.ipynb) for applying diHMM-cpp for H3K4me3 in hESC H1 Cells.

## Visualization the results in Genome browser

![example](images/H1_H3K4me3_CDC14A.jpg)





