# Bayesian Psychometric Network Analysis of Depression and OCD Symptoms

## Overview
This repository provides a comprehensive, reproducible framework for conducting **Psychometric Network Analysis (PNA)** on comorbid psychiatric symptoms. Using clinical data from the **Rogers Dataset**, this project implements and evaluates two distinct statistical paradigms for network estimation:
1. **Frequentist Regularized Networks** (via `EBICglasso`)
2. **Bayesian Graphical Models** (via Markov Chain Monte Carlo / Gibbs Sampling)

The primary empirical objective is to map out the interconnected symptom structures between **Depression** (QIDS-SR) and **Obsessive-Compulsive Disorder** (Y-BOCS-SR) in adult clinical patients, identifying central vulnerabilities and topological bridging pathways.

---

## Repository Structure
```text
├── scripts/
│   └── bayesian_psychometric_network.R  # Core analytical pipeline (Data prep, Frequentist, Bayes)
└── output/
    ├── frequentist_network.png          # EBICglasso regularized symptom network
    ├── centrality_plot.png             # Node Strength hierarchy analysis
    └── bayesian_network.png             # Bayesian posterior parameters network
```

---

## Methodological Blueprint

### 1. Psychometric Data Preparation
Clinical questionnaire data often exhibits artificial negative dependencies due to mutually exclusive filtering options. Using the adult sample ($N = 408$), this pipeline resolves structural redundancy within the QIDS-SR scale by collapsing sleep disturbances (items 1–4) and weight/appetite changes (items 6–9) into maximum-severity composite scores (`Sleep_Disturbance` and `Weight_Appetite`).

### 2. Frequentist Paradigm (`EBICglasso`)
* **Mathematical Core**: Partial correlation matrix regularized via the **Least Absolute Shrinkage and Selection Operator (LASSO)**.
* **Sparsity Optimization**: The Extended Bayesian Information Criterion (EBIC) with a tuning parameter ($\gamma = 0.5$) is applied to filter out statistical noise, driving trivial edges to absolute zero to prevent model overfitting.

### 3. Bayesian Paradigm (`easybgm` & `bgms`)
* **Mathematical Core**: A Bayesian Graphical Model designed for ordinal/Likert-scale data.
* **Sampling Engine**: Markov Chain Monte Carlo (MCMC) utilizing a Gibbs sampler to estimate the full joint posterior distribution.
* **Advantage over Frequentist**: Rather than applying a rigid, deterministic mathematical cutoff, the Bayesian framework calculates the **Posterior Inclusion Probability (PIP)** for every edge, explicitly quantifying parameter uncertainty and enabling the estimation of true **95% Credibility Intervals**.

---

## Key Clinical Findings

### Topological Properties & Subcommunity Separation
Both frameworks consistently demonstrate that symptoms cluster strongly into their expected nosological domains (Pastel Red for Depression; Pastel Blue for OCD). This reinforces the structural validity of the two distinct psychometric scales within a severe clinical cohort.

### Centrality Mechanisms (The Psychological Engine)
The frequentist Node Strength analysis revealed two primary drivers of the interconnected psychopathology system:
1. **`Interf_Comp` (Compulsion Interference)**: The functional impairment caused by physical behavioral rituals.
2. **`Mood` (Depressed Mood)**: The core affective component of depressive disorders.

Because `Mood` extends mild-to-moderate connections to nearly all surrounding symptoms (sleep issues, guilt, fatigue, suicidal ideation), its cumulative absolute edge weight drives it to the top of the centrality hierarchy, acting as a critical systemic engine.

### Comorbidity Bridging Pathways
The **Bayesian network matrix** captures highly nuanced, probabilistic edge distributions across the comorbidity gap. A highly resilient and theoretically vital bridge was mapped directly between **`DP3_Moo` (Mood)** and **`OCD10_CC` (Control of Compulsions)**, illustrating the explicit cognitive-affective pathway where emotional distress feeds into the loss of behavioral control in OCD.

---

## Requirements & Reproducibility
To replicate the entire pipeline, launch your R environment and execute:

```R
# Install required packages
install.packages(c("MPsychoR", "bootnet", "qgraph", "easybgm", "bgms", "dplyr"))

# Run the analytical script
source("scripts/bayesian_psychometric_network.R")
```
