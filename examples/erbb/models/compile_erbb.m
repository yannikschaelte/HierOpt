exdir = fileparts(which('erbb_standard.m'));

amiwrap('erbb_standard','erbb_standard_syms',exdir,1);
amiwrap('erbb_hierarchical_adjoint','erbb_hierarchical_adjoint_syms',exdir,1);