# HierOpt

A MatLab toolbox for the analytical hierarchical computation of scaling and noise parameters. 

Currently, it supports proportionality, offset and noise parameters for normally distributed errors.

## Requirements

The likelihood functions are defined so as to work directly with the AMICI ODE computation toolbox, and the test rely on AMICI and the PESTO parameter estimation toolbox. Both toolboxes are freely available at github.com/icb-dcm.
The HierOpt toolbox is however built in a modular way so that it is easy to define e.g. own likelihood functions, or extract functionality offered here and include it in an own optimization framework.

## Examples

In the examples, the hierarchical approach is used to perform parameter estimation on Systems Biology ODE models. To compute derivatives of the objective function, both the forward and adjoint approach can be used.