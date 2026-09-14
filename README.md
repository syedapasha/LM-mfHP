<div align="center">

# LM-mfHP

### Lagrange Multiplier Test for Event Covariates in Marked Hawkes Processes

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen?style=for-the-badge)](https://github.com/your-username/LM-mfHP/actions)
[![License](https://img.shields.io/github/license/your-username/LM-mfHP?style=for-the-badge)](LICENSE)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=for-the-badge)](https://github.com/your-username/LM-mfHP/pulls)
[![GitHub Stars](https://img.shields.io/github/stars/your-username/LM-mfHP?style=for-the-badge&color=yellow)](https://github.com/your-username/LM-mfHP/stargazers)

</div>

## Overview

> Event data, particularly in complex systems like financial markets, social networks, or seismic activity, often exhibit self-exciting behavior captured by Hawkes processes. However, effectively incorporating and statistically validating the influence of external factors (covariates) on these events, especially when marks are involved, presents a significant challenge. Existing methods are not general enough to robustly test for the significance of these covariates.

The `LM-mfHP` project is a Lagrange Multiplier test specifically designed for event covariates within a marked Hawkes process framework. This enables practitioners to determine the statistical significance of external influences on event occurrences and their associated marks, leading to hypothesis testing in complex event-driven systems.

## Features

Implements a robust Lagrange Multiplier test to assess the significance of event covariates in marked Hawkes processes, ensuring reliable statistical inferences.

Utilizes an Expectation-Maximization (EM) algorithm for efficient and accurate estimation of Hawkes process parameters, even with complex marked data.

Generate insightful plots for marked point processes and perform residual analysis to validate model fit and detect potential misspecifications.

Evaluate model efficacy through Receiver Operating Characteristic (ROC) curves, providing a clear measure of predictive power and classification accuracy.

Includes utilities for generating synthetic marked Hawkes process data, facilitating testing, simulation, and understanding of the model's behavior.

## Architecture

This project is developed in R, leveraging its statistical and data manipulation capabilities to implement the Lagrange Multiplier test and associated utilities.

### Directory Structure

```
LM-mfHP/
├── 📄 data_gen.r           # Script for generating synthetic marked Hawkes process data.
├── 📄 em.r                 # Implementation of the Expectation-Maximization (EM) algorithm for parameter estimation.
├── 📄 gen_mpp_plots.r      # Utility for generating plots specific to marked point processes.
├── 📄 main.r               # Main script to run the Lagrange Multiplier test and analysis.
├── 📄 main_mc.r            # Main script for Monte Carlo simulations.
├── 📄 plots.r              # General plotting utilities for results and diagnostics.
├── 📄 res2.RData           # Saved results data (e.g., from a simulation run).
├── 📄 resMM2.RData         # Another set of saved results data.
├── 📄 residual.r           # Script for performing and analyzing model residuals.
└── 📄 roc.r                # Script for generating Receiver Operating Characteristic (ROC) curves.
```

## Operational Setup

### Prerequisites

To run this project, you will need:

*   **R**: Version 4.0 or higher.
    *   Download and install R from [CRAN](https://cran.r-project.org/).
*   **RStudio (Recommended)**: An integrated development environment for R.
    *   Download and install RStudio Desktop from [RStudio](https://www.rstudio.com/products/rstudio/download/).

### Installation

Follow these steps to get `LM-mfHP` up and running on your local machine:

1.  **Clone the Repository**:
    Open your terminal or RStudio console and clone the GitHub repository:

    ```bash
    git clone https://github.com/your-username/LM-mfHP.git
    cd LM-mfHP
    ```

2.  **Install R Packages**:
    The project relies on several standard R packages. You can install them by running the following commands in your R console:

    ```R
    # Install necessary packages (example - adjust based on actual dependencies)
    install.packages(c("data.table", "ggplot2", "pROC", "stats", "utils"))
    # Add any other specific packages required by data_gen.r, em.r, etc.
    # For instance, if a specific Hawkes process package is used, it should be listed here.
    ```

3.  **Run the Main Script**:
    Once the packages are installed, you can execute the main analysis script:

    ```R
    source("main.r")
    # Or for Monte Carlo simulations:
    # source("main_mc.r")
    ```

    This will run the Lagrange Multiplier test, perform estimations, and generate relevant outputs and plots as defined within the scripts.
