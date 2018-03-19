exdir = fileparts(which('jakstat_hierarchical_forward.m'));

amiwrap('jakstat_standard','jakstat_standard_syms',exdir,1);
amiwrap('jakstat_hierarchical_forward','jakstat_hierarchical_forward_syms',exdir,1);
amiwrap('jakstat_hierarchical_forward_offsets','jakstat_hierarchical_forward_offsets_syms',exdir,1);
amiwrap('jakstat_hierarchical_adjoint','jakstat_hierarchical_adjoint_syms',exdir,1);
amiwrap('jakstat_hierarchical_adjoint_offsets','jakstat_hierarchical_adjoint_offsets_syms',exdir,1);