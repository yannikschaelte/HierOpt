[exdir,~,~]=fileparts(which('rafmekerk_standard_syms.m'));

amiwrap('rafmekerk_standard','rafmekerk_standard_syms',exdir)
amiwrap('rafmekerk_hierarchical','rafmekerk_hierarchical_syms',exdir);
amiwrap('rafmekerk_hierarchical_adjoint','rafmekerk_hierarchical_adjoint_syms',exdir);
amiwrap('rafmekerk_hierarchical_noreps','rafmekerk_hierarchical_noreps_syms',exdir);
amiwrap('rafmekerk_hierarchical_noreps_adjoint','rafmekerk_hierarchical_noreps_adjoint_syms',exdir);