%% compilation JakStat for standard

[exdir,~,~]=fileparts(which('JakStat_syms.m'));
% compile the model
amiwrap('JakStat','JakStat_syms',exdir)

%% compilation JakStat hierarchical

[exdir,~,~]=fileparts(which('JakStat_hierarchical_syms.m'));
% compile the model
amiwrap('JakStat_hierarchical','JakStat_hierarchical_syms',exdir)

amiwrap('JakStat_hierarchical_adjoint','JakStat_hierarchical_adjoint_syms',exdir)