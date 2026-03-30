# Prostate-Cancer-Analysis Project
This project improves prostate cancer prediction by combining PSA with factors like age, ethnicity, and deprivation. It compares logistic regression and machine learning models to assess accuracy and fairness in risk prediction across different populations.

![Status](https://img.shields.io/badge/status-active-brightgreen)
![Language](https://img.shields.io/badge/language-R%20%7C%20Python-blue)
![License](https://img.shields.io/badge/license-Academic-lightgrey)

## Description

This project improves prostate cancer prediction by combining PSA with factors like age, ethnicity, and deprivation. It compares logistic regression and machine learning models to assess accuracy and fairness in risk prediction across different populations.

## Overview

This project focuses on improving prostate cancer prediction using prostate-specific antigen (PSA) levels alongside other important risk factors. Current screening methods often rely on a single PSA threshold, which may not be accurate for all population groups. This study explores whether combining multiple variables can improve prediction accuracy and fairness. The analysis uses both statistical and machine learning approaches to better understand how different factors contribute to prostate cancer risk.

## Objectives

The main objective of this project is to evaluate the role of PSA within a multivariable prediction model. It aims to assess whether adding factors such as age, ethnicity, and socioeconomic status improves prediction performance. Another objective is to compare traditional statistical models, such as logistic regression, with machine learning methods like random forests. The project also examines whether these models provide more accurate and fair predictions across different population groups.

## Data Structure

The dataset includes patient-level information related to prostate cancer risk. Key variables include PSA levels, age, ethnicity, and measures of socioeconomic status such as deprivation. Additional clinical or lifestyle variables may also be included where available. The data is organised in a tabular format, where each row represents an individual and each column represents a variable. Date variables are standardised, and missing data is handled carefully to ensure consistency and reliability.

## Pseudocode

The analysis begins by importing the dataset and required libraries. Data cleaning is then carried out, including formatting dates, handling missing values, and checking variable types. Exploratory data analysis is performed to understand the distribution of variables and identify patterns. A baseline logistic regression model is fitted using PSA as the main predictor, followed by a multivariable model including additional risk factors. A machine learning model, such as a random forest, is then trained for comparison. Model performance is evaluated using metrics such as AUC, and results are compared to determine whether more complex models improve prediction accuracy and fairness.

## Project Structure

```
project-folder/
│── data/              # Raw and processed datasets (if allowed)
│── scripts/           # Analysis scripts (R or Python)
│── outputs/           # Figures, tables, results
│── report/            # Final report (PDF/Word)
│── README.md          # Project documentation
```

## How to Use

Clone the repository, open the script files in your preferred environment (e.g. RStudio), and run the analysis. Ensure all required packages are installed before execution.

## Author

Olaoluwa Adebayo
Exeter University
Health Data Analytics/ Module Name




