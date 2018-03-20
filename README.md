# HierOpt

A MatLab toolbox for the analytical hierarchical computation of scaling and noise parameters. 

Currently, it supports proportionality, offset and noise parameters for normally distributed errors.

## Overview

The toolbox is built in a modular way making it easy to extract part of the functionality and include it in an own framework:
In the hieropt_b, hieropt_c and hieropt_noise functions, optimal values for offset, proportionaltiy and noise parameters are computed. The data for these computations are extracted and organized in hieropt_scalings. More, hieropt_nllh_adjoint and hieropt_nllh_forward are likelihood functions which enable using HierOpt in combination with the ODE computation toolbox AMICI, where the adjoint and forward approaches can be used to compute sensitivities, which are often useful for optimization of the likelihood function.

## Requirements

The likelihood functions are defined so as to work directly with AMICI, and the examples use AMICI and the parameter estimation toolbox PESTO. Both toolboxes are freely available at [github.com/icb-dcm].
As mentioned before, due to the modularity of the HierOpt toolbox, it is however easy to replace the likelihood functions by own likelihoods, so that this toolbox has no real requirements and is generally applicable.

## Examples

In the examples, the hierarchical approach is used to perform parameter estimation on Systems Biology ODE models. To compute derivatives of the objective function, both the forward and adjoint approach can be used.