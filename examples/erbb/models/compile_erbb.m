exdir = fileparts(which('erbb_standard_syms.m'));

o2flag = 0;

amiwrap('erbb_standard','erbb_standard_syms',exdir,o2flag);
amiwrap('erbb_hierarchical_forward','erbb_hierarchical_forward_syms',exdir,o2flag);
amiwrap('erbb_hierarchical_adjoint','erbb_hierarchical_adjoint_syms',exdir,o2flag);