[exdir,~,~]=fileparts(which('rafmekerk_standard_syms.m'));

o2flag = 0;

amiwrap('rafmekerk_standard','rafmekerk_standard_syms',exdir,o2flag)
amiwrap('rafmekerk_hierarchical_forward','rafmekerk_hierarchical_forward_syms',exdir,o2flag);
amiwrap('rafmekerk_hierarchical_adjoint','rafmekerk_hierarchical_adjoint_syms',exdir,o2flag);
amiwrap('rafmekerk_hierarchical_noreps_forward','rafmekerk_hierarchical_noreps_forward_syms',exdir,o2flag);
amiwrap('rafmekerk_hierarchical_noreps_adjoint','rafmekerk_hierarchical_noreps_adjoint_syms',exdir,o2flag);

amiwrap('rafmekerk_laplace_standard','rafmekerk_standard_syms',exdir,o2flag)
amiwrap('rafmekerk_laplace_hierarchical_noreps_forward','rafmekerk_laplace_hierarchical_noreps_forward_syms',exdir,o2flag);
amiwrap('rafmekerk_laplace_hierarchical_noreps_adjoint','rafmekerk_laplace_hierarchical_noreps_adjoint_syms',exdir,o2flag);