%% compilation RafMekErk for standard
try 
    rmdir('/Users/carolinloos/PhD/Software/AMICIGit/models/RafMekErk','s')
catch
end
[exdir,~,~]=fileparts(which('RafMekErk_syms.m'));
% compile the model
amiwrap('RafMekErk','RafMekErk_syms',exdir)

%% compilation RafMekErk for hierarchical
try 
    rmdir('/Users/carolinloos/PhD/Software/AMICIGit/models/RafMekErk_hierarchical','s')
catch
end
[exdir,~,~]=fileparts(which('RafMekErk_hierarchical_syms.m'));
% compile the model
amiwrap('RafMekErk_hierarchical','RafMekErk_hierarchical_syms',exdir)

%% compilation RafMekErk for hierarchical with adjoints

[exdir,~,~]=fileparts(which('RafMekErk_hierarchical_adjoint_syms.m'));
amiwrap('RafMekErk_hierarchical_adjoint','RafMekErk_hierarchical_adjoint_syms',exdir);

[exdir,~,~]=fileparts(which('RafMekErk_hierarchical_adjoint2_syms.m'));
amiwrap('RafMekErk_hierarchical_adjoint2','RafMekErk_hierarchical_adjoint2_syms',exdir);

[exdir,~,~]=fileparts(which('RafMekErk_hierarchical_adjoint_reps_syms.m'));
amiwrap('RafMekErk_hierarchical_adjoint_reps','RafMekErk_hierarchical_adjoint_reps_syms',exdir,1);