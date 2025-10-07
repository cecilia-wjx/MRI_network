# Multi-Organ Network Analysis

This repository contains code used to generate the results of the article titled "Multi-Organ Network of Cardiometabolic Disease-Depression Multimorbidity Revealed by Imaging Phenotypic and Genetic Analyses".

## Analysis Pipeline

The analysis workflow follows the methodology described in our manuscript:

### 1. Phenotypic Association Analysis
- **`traits_disease_association.R`** - Identifies imaging-derived phenotypes from abdomen, heart, and brain significantly associated with CMDs-depression multimorbidity using multivariable logistic regression.

### 2. Cross-Organ Correlation Analysis
- **`cross_organ_association.R`** - Quantifies pairwise correlations among multimorbidity-associated imaging traits across different organ systems using multivariable linear regression, constructing abdomen-heart-brain cliques.

### 3. Genetic Architecture Analysis

#### (1) Genetic Correlation
- **`LDSC_heart_brain.sh`** - Estimates genome-wide genetic correlations between imaging traits using linkage disequilibrium score regression (LDSC).
- **`lava.R`** - Performs local genetic correlation analysis to identify regional pleiotropic effects across abdomen-heart-brain cliques.
- **`GPA.R`** - Conducts genetic analysis incorporating pleiotropy and annotation (GPA) to estimate variant-level pleiotropy and genetic overlap between imaging traits.

#### (2) Shared Genomic Loci Identification
- **`cpassoc.R`** - Identifies shared genomic loci across abdomen-heart-brain cliques using cross-phenotypic association analysis (CPASSOC).

#### (3) Transcriptome-Wide Association
- **`TWAS.sh`** - Conducts transcriptome-wide association studies to identify genes expressed in organ-related tissues that are shared across abdomen-heart-brain cliques.

### 4. Multimorbidity Risk Prediction
- **`Prediction_Model.R`** - Develops and evaluates Cox proportional hazards models incorporating genetically-informed multi-organ imaging traits to predict multimorbidity development and progression.
