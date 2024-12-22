# Survival Analysis of Treatment Efficacy in Primary Biliary Cirrhosis (PBC)

## Overview

This project analyzes the survival outcomes of patients with Primary Biliary Cirrhosis (PBC) based on treatment with Cyclosporin A (CyA) and other clinical covariates. 
The analysis utilizes survival analysis techniques such as Kaplan-Meier estimators, Cox Proportional Hazards Models, and Accelerated Failure Time (AFT) models to evaluate treatment efficacy and identify key prognostic factors.

## Key Features

- **Dataset**: The PBC3 dataset from a multicenter randomized clinical trial conducted between 1983-1987 across six European hospitals.
- **Primary Outcome**: Time to failure of medical treatment (encompassing death or liver transplantation).
- **Analysis Methods**:
  - Descriptive statistics and visualizations.
  - Kaplan-Meier survival analysis.
  - Regression modeling using Cox Proportional Hazards and AFT models.
  - Proportional hazards assumption checks.

## Results

- No statistically significant improvement in survival for patients treated with CyA (p = 0.7813).
- Significant predictors of survival:
  - **Sex**: Males had higher hazard ratios.
  - **Bilirubin**: Higher levels associated with shorter survival.
  - **Albumin**: Higher levels associated with longer survival.
  - **Disease Stage**: Advanced stages linked to poorer outcomes.
- Kaplan-Meier analysis revealed significant differences based on sex, disease stage, and gastrointestinal bleeding history.

## Highlights

- Comprehensive diagnostic checks confirmed the robustness of the Cox and Weibull AFT models.
- Key insights on time-dependent effects of disease stages and interactions between covariates.

## Repository Structure

- `data/`: Includes the PBC3 dataset (access restricted).
- `code/`: SAS code files for data analysis and visualization.
- `results/`: Outputs from analyses, including Kaplan-Meier plots, regression results, and diagnostics.
- `docs/`: Documentation and references.
