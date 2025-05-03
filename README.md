# Metabolic Regulatory Network
This repository contains code and resources for modeling and analyzing metabolic regulatory networks.

## Contents
- **bioutils/**: Utilities for biological computations.
- **MRN/**: Core metabolic regulatory network functions and models.
- Various MATLAB scripts for different analysis workflows:
  - `MRN_main.m`
  - `S1_MRN.m` to `S8_MRN.m`
  - and others.

## Usage
To run the models, open the desired `.m` file in MATLAB and execute.

## Requirements
- MATLAB R2021b or newer.
- Required MATLAB toolboxes:
  - Symbolic Math Toolbox
  - Optimization Toolbox

## License
Add the license information here.






# Metabolic Regulatory Network

This repository contains code and resources for modeling and analyzing **Metabolic Regulatory Networks (MRNs)**, with advanced capabilities for addressing biological complexity using **continuous models** and **sparse system identification**.

## Overview

Traditional hybrid models using Boolean logic become impractical as the number of regulatory proteins increases exponentially. This repository introduces a **continuous modeling framework** that replaces Boolean rules with **Hill functions**, enabling more realistic and computationally efficient simulations of metabolic-regulatory networks.

Additionally, we integrate a **Sparse Identification of Nonlinear Dynamical Systems (SINDy)** framework, which infers compact, interpretable models directly from simulation or empirical data. This enables automatic discovery of governing equations and reduces the need for manually specified parameters.

## Highlights

- 📈 **Continuous Model with Hill Functions:** Simplifies regulatory interactions and reduces computational complexity.
- 📉 **Sparse System Identification (SINDy):** Efficiently infers the governing equations of MRNs using rational functions and polynomial terms.
- ⚡ **Reduced Model Complexity:** Requires fewer parameters than hybrid models while capturing essential nonlinear dynamics.
- 🔬 **Biologically Realistic Behavior:** Captures graded responses of regulatory proteins and metabolites.

## Contents

- **bioutils/**: Utilities for biological computations.
- **MRN/**: Core metabolic regulatory network functions and models.
- Various MATLAB scripts for different analysis workflows:
  - `MRN_main.m`: Main entry point for simulations.
  - `S1_MRN.m` to `S8_MRN.m`: Various experimental workflows and case studies.
- Python (optional): Sparse system identification workflow implemented externally (not included here).

## Usage

To run the models:

1. Open the desired `.m` file in **MATLAB**.
2. Execute the script to simulate metabolic and regulatory dynamics.
3. Use the SINDy-based identification workflow (optional, in Python) to infer reduced models from the simulation results.

## Requirements

- MATLAB R2021b or newer.
- Required MATLAB Toolboxes:
  - Symbolic Math Toolbox
  - Optimization Toolbox

Optional (for Sparse Identification / Python version):

- Python 3.x
  - numpy
  - scipy
  - scikit-learn (or other sparse regression libraries)

## Features and Research Contributions

- **Continuous Model**:
  - Replaces binary logic with smooth regulatory dynamics using Hill functions.
  - Captures nuanced biological processes including gene expression regulation and metabolite consumption.

- **Sparse System Identification (SINDy)**:
  - Discovers governing equations from time-series simulation or experimental data.
  - Supports rational functions and implicit ODE formats.
  - Reduces model size while retaining biological interpretability.
  
- **Validation and Comparison**:
  - Benchmark against hybrid models with varying regulatory weight parameters (`βp` values).
  - Analyze simulation results to evaluate model accuracy and robustness.

## License

Add license information here.

## Citation

If you use this model or approach, please cite:

> Sparse Identification and Mathematical Framework for Analyzing Metabolic-Regulatory Networks  
> Neveen Ali Eshtewy, Ali Forootani, Shumaila Noreen, Mohammad Khosravi, 2025.

---

For more detailed explanations, theory, and results, refer to the full article and data:
- [Zenodo Archive](https://zenodo.org/records/14540008)
- [GitHub Repository](https://github.com/Ali-Forootani/Metabolic_regulatory_network)


