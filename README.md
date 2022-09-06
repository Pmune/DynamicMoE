# DynamicMoE

## Description

Dynamic Mixture of Experts models are flexible models for modeling the predictive density of a response variable as a mixture of experts model with covariate-dependent mixture components (experts) and mixture weights (gates). The parameters in both the components and the mixture weights are allowed to flexibly change over time.

The component models or experts can be any well-behaved density. By well-behaved density, we mean a density function that is twice deferentiable with respect to the model's parameter and that is identifiable.

This package implements an efficient particle filtering inference methodology for sampling from the posterior sequentially over time and for making online predictions. The algorithm is designed to handle static or dynamic mixture of experts models with high-dimensional parameter in a unified way.

The current implementation includes only Poisson component models and a single-compmonent generalized Poisson model. Training dynamic mixture of Poisson models is done through the function `dmoe`. The single-compmonent generalized Poisson model training is done through the function `dgenpois`.

## Installation

The package can be downloaded as:
```
library("devtools")
devtools::install_github(repo = "Pmune/DynamicMoE/dmoe")
```

## File organization

- The R package is saved in the folder `dmoe`. All codes are saved in the subfolder `R`.
- The folder `simulations` includes some examples.
## Examples

Examples are provided through simulation experiments, where data are generated from:

- M1: Static Poisson regression data generating process model.
- M2: Dynamic Poisson regression data generating process model.
- M3: Dynamic mixture of Poisson regression experts data generating process model.
- M4: Static mixture of Poisson autoregressive experts data generating process model.
- M5: Dynamic mixture of Poisson autoregressive experts data generating process model.

The aim of the above experiments is to check if the Sequential Monte Carlo (SMC) inference method implemented in the package `dmoe` can recognize the true underlying data generating process.
To do this,  $50$ datasets are generated from each data generating process and dynamic mixture of Poisson experts models with different parameter settings are trained. Among the trained models, one model is selected using log predictive scores as the model selection criterion. The best model is the one which has the highest log predictive score.

In the experiments M4 and M5, the SMC inference is compared with the MCMC method which is implemented in Matlab. The aim here is to see whether training dynamic models has better value over training static models.

## Instructions to run the simulation experiments

The files in the `simulations` folder implement simulation experiments. The simulations folder includes:

- The `R` folder. All source codes implementing the simulation experiments are saved in this folder.
- The `data` folder. Two csv files simdata_static.csv and simdata_dyn.csv, which contains $50$ datasets generated from M4 and M5 respectively, are saved in this folder.
- The `results` folder. All results are saved in this folder. In this folder includes also files containing the log predictive scores (LPS), computed from the MCMC sampler, required to complete the SMC and MCMC comparison.
- The `plots` folder. It contains plots and the code generating the plots.

The `R` folder includes the file `run_all_simulation_experiments.R` which automates all experiments. The code in this file does:

- Install the dmoe package
- Install all libraries required to run the simulation experiments and to generate the plots.
- Run each of the M1-M4 experiments and the generate the plots automatically.

To run the code in `run_all_simulation_experiments.R`:

- If you use Rstudio, double click on the simulations.Rproject. It will open up an Rstudio window with all files.
- Otherwise, set the working directory to the `simulations` folder and run the code in `run_all_simulation_experiments.R`.
