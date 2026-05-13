# Multi-scale Analysis of Taxonomic Decay and Functional Reshuffling

This repository contains the updated R script (`análises10.05.R`) used to investigate the decoupling between taxonomic and functional dimensions in dung beetle communities (Scarabaeinae) within the Atlantic Forest.

## 📋 Overview

The script implements a robust ecological workflow designed to test the **"Landscape Sovereignty"** and **"Neighborhood Determinism"** hypotheses. It focuses on how land-use intensification and soil properties drive the functional reshuffling ($r_s$) of communities at different spatial scales.

## 🚀 Workflow Description

The analysis is structured into six main logical blocks:

### 1. Environment and Reproducibility
- Loads essential ecological and statistical packages (`tidyverse`, `vegan`, `sads`, `MuMIn`, `logistf`, etc.).
- Uses the `here` package for robust relative path management.
- Sets a global seed to ensure reproducibility of stochastic elements (like jittering in TADs).

### 2. Data Synchronization (L, Q, R Matrices)
- Imports abundance, trait, and environmental data.
- **Key Step:** A custom standardization function cleans and synchronizes species names across all matrices, ensuring data integrity before downstream analyses.

### 3. Assembly Rule Diagnosis (SAD vs. TAD)
- Fits **Species Abundance Distributions (SAD)** and **Trait Abundance Distributions (TAD)**.
- Compares multiple models (Log-series, Log-normal, Power-law) via $\Delta$AIC to detect shifts from deterministic to stochastic signatures.

### 4. Scale of Effect (SoE) "Tournament"
- Implements a vectorized selection process to identify the optimal spatial radius for each environmental predictor.
- Tests multiple buffers (from **210m to 990m**) to find the scale where the biological response (taxonomic or functional) is strongest.

### 5. Statistical Modeling (LRT & Firth Regression)
- Conducts **Likelihood Ratio Tests (LRT)** to compare the relative importance of landscape vs. microclimatical drivers.
- **Firth’s Penalized Likelihood (`logistf`):** Applied to handle rare events and data separation issues, ensuring unbiased estimates even with low sample sizes in matrix habitats.
- Focuses on the **Functional Reshuffling ($r_s$)** of biomass and mechanical traits (protibia and metatibia).

### 6. High-Resolution Data Visualization
- Generates 600 DPI publication-quality panels using `ggplot2` and `patchwork`.
- Includes internal legends and theme customization to highlight the divergence between forest and matrix responses.

## 🛠 Requirements

To run this script, you need R installed with the following libraries:

```R
install.packages(c("tidyverse", "sads", "MuMIn", "openxlsx", "logistf", "vegan", "patchwork", "here"))
