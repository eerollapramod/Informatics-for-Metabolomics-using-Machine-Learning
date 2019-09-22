# Classification and Regression models using Machine Learning (Informatics for Metabolomics)

This project contains three seperate R scripts; 
  1. Data Preperation and Exploratory Analysis
  2. Classification models
  3. Regression models


## Requirements

You will need the `R 3.6.0` or later for some packages to run.

#### Dataset
Four separate CSV files: `Bacterial_Counts.csv, Enose_data.csv, HPLC_data.csv and Sensory_score.csv` are attached within the repository. These data were collected during an experiment in which fresh mince beef fillets were stored under different temperatures (0°C, 5°C, 10°C, 15°C, 20°C). Sampling was taking place at regular time intervals and the samples were subjected to enose profiling, solvent extraction and HPLC analysis and plate counting to account for Total Viable Counts (TVC) and Pseudomonas. In addition a taste panel was scoring the samples based on visual appearance and odour.

Sample ids (e.g. 0F4) are unique and can be used to match rows across the four files. The first number corresponds to the temperature, “F” stands for fillet, and the last number represents the sampling time point (hours). “F1a” and “F1b” are the two replicates sampled at baseline (i.e. on arrival to the lab).

The column names in the HPLC_data.csv represent the major peaks present in the samples which are denoted by their retention time (min). In the Sensory _score.csv the column ‘sensory’ is a categorical variable representing the sensory score attributed by the taste panel. In this case, sensory class 1 denotes fresh, class 2 denotes semi-fresh (just acceptable) and 3 denotes spoilage (1 = fresh, 2=semi-fresh, 3=spoiled).


## 1. Data Preparation and Exploratory analysis

In this section (`DataAnalysis.R` script), the data has been sorted by matching the identical sample ids and the following steps were performed: 
  * PCA plots for Enose and HPLC data
  * Generates the histograms for Total Viable Counts (TVC) during fresh, semi fresh and spoiled stages
  * Generates the histograms for Pseudomonas during fresh, semi fresh and spoiled stages
  * Plots the bacterial counts (TVC) against time for each temperature.
  * Plots the Pseudomonas count against time for each temperature.
  * Produces scatterplots for different combinations of sensors (for enose) or metabolites (for HPLC) and bacterial counts.
  * Produces boxplots of bacterial counts (TVC) against the sensory score.
  * Produces boxplots of Pseudomonas against the sensory score.


## 2. Classifiscation Models

This part (`classification.R` script) produces classification models for sensory scores using Enose data and HPLC data. Decision Trees, Standard Vector Machines, and k-Nearst Neighbours were used for clasiification. FOr each model, it generates:
* Generates evaluation metrics including ROC curves and confusion matrix (Uses sensory score=3 as the positive class, when producing the confusion matrix).
* Measures the success rate of each model over multiple iterations and presents the results in a plot for comparison purposes where necessary.
* Create a plot showing the importance of each variable for the prediction (if applicable).

## 3. Regression models

In part three, `regression.R` script generates threeregression models (Linear Model(LM), kNN, Random Forests) for TVC & Pseudomonas Counts by Enose data, and TVC & Pseudomonas Counts by HPLC data. For each model:
* It calculates the RMSE for each model.
* Produces plots comparing the observed against the predicted values of Total Viable Counts
(TVC) and Pseudomonas obtained by different regression methods and different analytical
platforms (In the plots, includes a fitted line and 95% confidence intervals (CI)).









