exdir = fileparts(which('jakstat_hierarchical_forward_syms.m'));

o2flag = 0;

amiwrap('jakstat_standard','jakstat_standard_syms',exdir,o2flag);
amiwrap('jakstat_hierarchical_forward','jakstat_hierarchical_forward_syms',exdir,o2flag);
amiwrap('jakstat_hierarchical_forward_offsets','jakstat_hierarchical_forward_offsets_syms',exdir,o2flag);
amiwrap('jakstat_hierarchical_adjoint','jakstat_hierarchical_adjoint_syms',exdir,o2flag);
amiwrap('jakstat_hierarchical_adjoint_offsets','jakstat_hierarchical_adjoint_offsets_syms',exdir,o2flag);

amiwrap('jakstat_laplace_standard','jakstat_laplace_standard_syms',exdir,o2flag);
amiwrap('jakstat_laplace_hierarchical_adjoint','jakstat_laplace_hierarchical_adjoint_syms',exdir,o2flag);
amiwrap('jakstat_laplace_hierarchical_adjoint_offsets','jakstat_laplace_hierarchical_adjoint_offsets_syms',exdir,o2flag);
