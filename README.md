# The Connecting Thread Theory (TBP): Computational Evaluation Framework

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![SymPy](https://img.shields.io/badge/SymPy-Fractional_Precision-green.svg)](https://www.sympy.org/)

## Overview
This repository provides the open-source computational evaluation script for the Connecting Thread Theory (TBP). The TBP framework proposes a 6-dimensional topological matrix modeled to evaluate fundamental physical constants. Within the limits of this geometric model, parameters such as particle mass hierarchies and cosmological energy budgets are analytically derived, structurally minimizing reliance on empirical curve-fitting. This script utilizes symbolic computation to verify the internal dimensional consistency and the mathematically singularity-free nature of the proposed framework.

## Explicit Assumptions and Limitations
The computations provided in this repository are strictly theoretical. The primary assumption of this framework is the formalization of a quantized spatial capacity constrained by an absolute thermodynamic ceiling ($R_{crit} \approx 215.11$ GeV). The script provides evaluations up to first-order geometric friction and second-order kinematic backscattering limits. Extrapolations beyond these evaluated parameters, particularly concerning higher-order quantum field fluctuations, remain speculative and are subject to future methodological development.

## Computational Stages
The codebase executes a linear, 6-stage algorithmic verification of the TBP mathematical architecture:

1. **Fundamental Parameters & Electroweak Sector**: Evaluates matrix dimensionality and formally derives the theoretical Higgs scalar equilibrium ($\sim 125.22$ GeV).
2. **QCD & Cosmology**: Models the proton equilibrium ($\sim 938.26$ MeV) and calculates a theoretical base for the Dark Sector (evaluating Dark Energy at $\approx 68.67\%$).
3. **Fermion Spectrum**: Derives lepton and quark mass limits using deterministic topological friction operators.
4. **Boundary Limits & Unitarity**: Evaluates high-energy asymptotic bounce limits and the fundamental $SU(2)$ isospin asymmetry ratio ($17/18$).
5. **Extreme Conditions**: Models the potential isolation of structural topological defects when energy inputs exceed $R_{crit}$.
6. **Higher-Order Perturbations**: Computes 2-loop kinematic backscattering and the localized terrestrial observer bias ($\Sigma \Phi_{local}$).

## Reproducibility and Usage
To ensure absolute computational reproducibility, the theoretical derivations can be executed and verified using the standard Python environment. The script relies exclusively on the core `math` library and `SymPy` to maintain exact fractional precision.

**Prerequisites:**
```bash
pip install sympy
```

**Execution:**
```bash
python tbp_computational_evaluation.py
```

## Data Availability
The empirical values of fundamental physical constants utilized for comparative analysis in this theoretical investigation are referenced from the Particle Data Group (PDG) and CODATA. This repository contains only the theoretical derivation algorithms; no new primary empirical datasets were generated.

## Citation and License
This codebase is an integral supplement to the theoretical formalization of the Connecting Thread Theory. It is released under the **Creative Commons Attribution-NonCommercial (CC BY-NC 4.0)** license to support open-science methodological transparency.

**Author:** Cefiyana  
**ORCID:** [https://orcid.org/0009-0008-4324-9515](https://orcid.org/0009-0008-4324-9515)
