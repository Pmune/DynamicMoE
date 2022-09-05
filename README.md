# DynamicMoE

## Description

Dynamic Mixture of Experts models are flexible models for modeling the predictive density of a response variable as a mixture of experts model with covariate-dependent mixture components (experts) and mixture weights (gates). The parameters in both the components and the mixture weights are allowed to flexibly change over time.

The component models or experts can be any well-behaved density. By well-behaved density, we mean a density function that is twice deferentiable with respect to the model's parameter and that is identifiable.

This package implements an efficient particle filtering inference methodology for sampling from the posterior sequentially over time and for making online predictions. The algorithm is designed to handle static or dynamic mixture of experts models with high-dimensional parameter in a unified way.

The current implementation includes only Poisson component models and a single-compmonent generalized Poisson model. Training dynamic mixture of Poisson models is done through the function `dmoe`. The single-compmonent generalized Poisson model training is done through the function `dgenpois`.

## Installation

The package can be downloaded as a library. To do this:
- Install the `devtools` library.

```
install.library("devtools") # if not installed
```
- Install the `dmoe` package from the Github repo `Pmune/DynamicMoE`

```
library("devtools")
devtools::install_github(repo = "Pmune/DynamicMoE/dmoe")
```
