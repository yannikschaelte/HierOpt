# HierOpt

A MatLab toolbox for the analytical hierarchical computation of scaling and noise parameters. 

Currently, it supports proportionality, offset and noise parameters for normally distributed errors.

## Overview and Usage

HierOpt is intended to be used when you want to optimize a parametrized likelihood function, some of whose parameters are scaling factors that do not influence your simulations. If the structure of your likelihood corresponds to one supported by this toolbox, you can in your likelihood define an inner loop in which you optimize the scaling parameters. That means, you first perform unscaled simulations, then compute optimal scalings, and then use these to compute your full likelihood and its parameters, based on the scaled simulations.

The toolbox is built in a modular way making it easy to extract part of the functionality and include it in an own framework:

In the hieropt_b, hieropt_c and hieropt_noise functions, optimal values for offset, proportionaltiy and noise parameters are computed. In the header of this function, you can also see the options struct that you have to create, specifying how the scalings shall be computed. The data for these calculations are extracted and organized in [hieropt_scalings](hieropt_scalings.m). More, [hieropt_nllh_adjoint](hieropt_nllh_adjoint.m) and [hieropt_nllh_forward](hieropt_nllh_forward.m) are likelihood functions which enable using HierOpt in combination with the ODE computation toolbox AMICI, where the adjoint and forward approaches can be used to compute sensitivities, which are often useful for optimization of the likelihood function. Here, you only need to pass a function handle to your simulation function as well as your data, then HierOpt will compute the full likelihood for you.

## Requirements

The likelihood functions are defined so as to work directly with AMICI, and the examples use AMICI and the parameter estimation toolbox PESTO. Both toolboxes are freely available at [github.com/icb-dcm](https://github.com/icb-dcm).
As mentioned before, due to the modularity of the HierOpt toolbox, it is however easy to replace the likelihood functions by own likelihoods, so that this toolbox has no real requirements and is generally applicable.

## Examples

In the examples, the hierarchical approach is used to perform parameter estimation on Systems Biology ODE models. The considered models are

* [jakstat_swameye](examples/jakstat_swameye): a model of the JakStat signaling pathway with time-resolved input ofthe drug EPO,
* [rafmekerk](examples/rafmekerk): a model of the RafMekErk signaling pathway,
* [erbb](examples/erbb):  a model of ErbB signaling pathways.

To compute derivatives of the objective function, both the forward and adjoint approach can be used in each of the models. The forward approach is sometimes faster in small-scale models, whereas for large-scale models the adjoint approach can be considerably more efficient.