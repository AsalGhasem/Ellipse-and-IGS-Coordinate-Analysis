# Ellipse and IGS Coordinate Analysis

## Overview
This repository contains MATLAB scripts and datasets for analyzing elliptical data and IGS (International GNSS Service) coordinates. The project is divided into two primary modules—`ellipse` and `igs`—both of which utilize Chi-squared statistical tests for data validation and analysis.

## Repository Structure

The project is organized under the main `pr1/` directory:

*   **`pr1/ellipse/`**
    *   `Data.mat`: MATLAB data workspace file containing the primary datasets for ellipse modeling.
    *   `chi2Test.m`: MATLAB function/script for performing Chi-squared tests on the ellipse data.
    *   `part1.m`: Execution script for the first part of the project.
    *   `part3.m`: Execution script for the third part of the project.

*   **`pr1/igs/`**
    *   `igs_coordinates.csv`: Raw CSV dataset containing the IGS coordinate data.
    *   `chi2Test.m`: MATLAB function/script for performing Chi-squared tests on the IGS coordinates.
    *   `part2.m`: Execution script for the second part of the project, focusing on IGS data.

## Prerequisites
*   **MATLAB**: Required to execute the `.m` scripts and load the `.mat` data files.

## Usage
1. Clone this repository to your local machine.
2. Open MATLAB and set your Current Folder to either `pr1/ellipse/` or `pr1/igs/` depending on the module you want to run.
3. Execute the scripts in their numbered order (`part1.m`, `part2.m`, `part3.m`) from the MATLAB command window or editor.
