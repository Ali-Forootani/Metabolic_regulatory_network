# Metabolic Regulatory Network

This repository contains code and resources for modeling and analyzing **Metabolic Regulatory Networks (MRNs)** using a continuous modeling approach combined with sparse identification techniques for discovering nonlinear dynamics.

## Overview

Traditional hybrid models using Boolean logic become increasingly complex and computationally expensive as the number of regulatory proteins grows. This repository introduces a **continuous mathematical model** which replaces Boolean logic with **Hill functions** to simplify regulatory dynamics and enable more biologically realistic simulations.

Furthermore, we integrate a **Sparse Identification and Mathematical Framework for Analyzing Metabolic-Regulatory Networks**, which allows for efficient discovery of governing equations directly from simulation or empirical data. This approach reduces the need for predefined model structures and identifies compact, interpretable models.

## Highlights

- 📈 **Continuous Model with Hill Functions:** Provides smooth and biologically relevant representations of regulatory interactions.
- 📉 **Sparse Identification and Mathematical Framework:** Efficiently discovers nonlinear dynamical models from simulation data.
- ⚡ **Reduced Model Complexity:** Simplifies formulation with fewer parameters than hybrid models.
- 🔬 **Biologically Realistic Behavior:** Captures graded regulatory responses and dynamic transitions in metabolic-regulatory systems.

## Contents

- **bioutils/**: Utilities for biological computations.
- **MRN/**: Core metabolic regulatory network models and functions.
- MATLAB analysis and simulation scripts:
  - `MRN_main.m`: Main entry point script for running models.
  - `S1_MRN.m` to `S8_MRN.m`: Experimental scenarios and analysis workflows.

## Usage

1. Open the desired `.m` file in **MATLAB**.
2. Execute the script to simulate the metabolic-regulatory network using the continuous model.
3. (Optional) For advanced analysis and model discovery, apply the **Sparse Identification and Mathematical Framework for Analyzing Metabolic-Regulatory Networks** using data from simulations.

## Requirements

- MATLAB R2021b or newer.
- MATLAB Toolboxes:
  - Symbolic Math Toolbox
  - Optimization Toolbox

(Optional Python for model discovery):

- Python 3.x
  - numpy
  - scipy
  - scikit-learn (or similar for sparse regression)

## Features and Research Contributions

- **Continuous Model for MRNs**:
  - Utilizes Hill functions to replace Boolean logic for regulatory proteins.
  - Applies mass action and Michaelis-Menten kinetics for metabolic and degradation processes.
  - Achieves smooth, realistic dynamics in regulatory interactions.

- **Sparse Identification and Mathematical Framework for Analyzing MRNs**:
  - Automatically infers governing equations directly from time series data.
  - Supports identification of nonlinear dynamics involving rational and polynomial terms.
  - Offers interpretable and reduced representations of system dynamics.

- **Model Validation**:
  - Compares continuous model results against hybrid model simulations.
  - Adjusts regulatory weight parameters (βp values) to analyze dynamic behavior under various conditions.

## License

Add license information here.

## Citation

If you use this repository or methodology, please cite:

> Sparse Identification and Mathematical Framework for Analyzing Metabolic-Regulatory Networks  
> Neveen Ali Eshtewy, Ali Forootani, Shumaila Noreen, Mohammad Khosravi, 2025.
> Email: Aliforootani@ieee.org

---

For further information, detailed theory, and data, refer to:

- [Zenodo Archive](https://zenodo.org/records/14540008)
- [GitHub Repository](https://github.com/Ali-Forootani/Metabolic_regulatory_network)


