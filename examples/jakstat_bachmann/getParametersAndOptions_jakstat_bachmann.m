function [parameters,options] = getParametersAndOptions_jakstat_bachmann(approach)

nStarts = 20;

options.MS = PestoOptions();
options.MS.n_starts = nStarts; % actually 500
options.MS.mode = 'text';
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
    'display','iter',...
    'GradObj','on',...
    'MaxIter',5000,...
    'TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxFunEvals',40000,...
    'PrecondBandWidth', inf);
options.MS.obj_type = 'negative log-posterior';

nPar = 112; % maximum number of parameters
% 27 dyn, 31 offsets, 43 proportionality, 11 sd
minPar = -3*ones(nPar,1);
maxPar = 3*ones(nPar,1);

maxPar(1)  = 4;  % CISEqc
maxPar(3)  = 12; % CISInh
maxPar(7)  = 4;  % EpoRActJAK2
maxPar(8)  = 6;  % EpoRCISInh
maxPar(10) = 9;  % JAK2ActEpo
maxPar(11) = 4;  % JAK2EpoRDeaSHP1
maxPar(20) = 4;

minPar(28:end) = -5; % offsets

rng(0);
par0 = bsxfun(@plus,minPar,bsxfun(@times,maxPar - minPar, lhsdesign(nStarts,nPar,'smooth','off')'));
   
switch approach
    case 'standard'
        nPar = 112;
        
    case 'hierarchical-adjoint'
        nPar = 58;
        
        sc.exp_groups.bc_idxs = {1};
        sc.exp_groups.noise_idxs = {1};
        
        sc.obs_groups.bc_idxs = {1,2,3};
        sc.obs_groups.b_mode = {'absolute','absolute','absolute'};
        sc.obs_groups.c_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.noise_idxs = {1,2,3};
        sc.obs_groups.noise_mode = {'multiple','multiple','multiple'};
        
        sc.distribution = 'normal';
        
        options.sc = sc;
        
    case 'hierarchical'
        nPar = 58;
        
        sc.exp_groups.bc_idxs = {1,2,3,[4,5],6,[7,8],[9,10],[11,12],[13,14],15:19,20:25,26:31,32:36};
        sc.exp_groups.noise_idxs = {1:36};
        
        for iyg = 1:20
            sc.obs_groups.bc_idxs{iyg} = iyg;
        end
        for iyg = [1:6,12:20]
            sc.obs_groups.b_mode{igy} = 'absolute';
            sc.obs_groups.c_mode{igy} = 'single';
        end
        for iyg = {7:11}
            sc.obs_groups.b_mode{igy} = 'absolute';
            sc.obs_groups.c_mode{igy} = 'absolute';
        end
        sc.obs_groups.noise_idxs = {[1,2],[3,19,20],4,[5,6],7,8,9,10,11,[12,13,14,15,16,17],18};
        for iyg = 1:11
            sc.obs_groups.noise_mode{igy} = 'single';
        end
        
        sc.distribution = 'normal';
        
        options.sc = sc;
        
    case 'hierarchical-offsets'
        nPar = 27;

        sc.exp_groups.bc_idxs = {1};
        sc.exp_groups.noise_idxs = {1};
        
        sc.obs_groups.bc_idxs = {1,2,3};
        sc.obs_groups.b_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.c_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.noise_idxs = {1,2,3};
        sc.obs_groups.noise_mode = {'multiple','multiple','multiple'};
        
        sc.distribution = 'normal';
        
        options.sc = sc;
        
    case 'hierarchical-adjoint-offsets'
        nPar = 27;
        
        
        sc.exp_groups.bc_idxs = {1};
        sc.exp_groups.noise_idxs = {1};
        
        sc.obs_groups.bc_idxs = {1,2,3};
        sc.obs_groups.b_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.c_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.noise_idxs = {1,2,3};
        sc.obs_groups.noise_mode = {'multiple','multiple','multiple'};
        
        sc.distribution = 'normal';
        
        options.sc = sc;
        
end

parameters.name = {...
    'CISEqc'
    'CISEqcOE'
    'CISInh'
    'CISRNADelay'
    'CISRNATurn'
    'CISTurn'
    'EpoRActJAK2'
    'EpoRCISInh'
    'EpoRCISRemove'
    'JAK2ActEpo'
    'JAK2EpoRDeaSHP1'
    'SHP1ActEpoR'
    'SHP1Dea'
    'SHP1ProOE'
    'SOCS3Eqc'
    'SOCS3EqcOE'
    'SOCS3Inh'
    'SOCS3RNADelay'
    'SOCS3RNATurn'
    'SOCS3Turn'
    'STAT5ActEpoR'
    'STAT5ActJAK2'
    'STAT5Exp'
    'STAT5Imp'
    'init_EpoRJAK2'
    'init_SHP1'
    'init_STAT5'
    'offset_CIS_actd'
    'offset_CIS_cisoe'
    'offset_CIS_long'
    'offset_CIS_shp1oe'
    'offset_CIS_socs3oe'
    'offset_SOCS3_cisoe'
    'offset_SOCS3_long'
    'offset_SOCS3_socs3oe'
    'offset_pEpoR_actd'
    'offset_pEpoR_cisoe'
    'offset_pEpoR_cisoe_pepor'
    'offset_pEpoR_dr30'
    'offset_pEpoR_dr7'
    'offset_pEpoR_fine'
    'offset_pEpoR_long'
    'offset_pEpoR_shp1oe'
    'offset_pEpoR_socs3oe'
    'offset_pJAK2_actd'
    'offset_pJAK2_cisoe'
    'offset_pJAK2_dr30'
    'offset_pJAK2_dr7'
    'offset_pJAK2_fine'
    'offset_pJAK2_long'
    'offset_pJAK2_shp1oe'
    'offset_pJAK2_socs3oe'
    'offset_pSTAT5_actd'
    'offset_pSTAT5_cisoe'
    'offset_pSTAT5_conc'
    'offset_pSTAT5_long'
    'offset_pSTAT5_shp1oe'
    'offset_pSTAT5_socs3oe'
    'scale1_CIS_dr90'
    'scale2_CIS_dr90'
    'scale_CISRNA_foldA'
    'scale_CISRNA_foldB'
    'scale_CISRNA_foldC'
    'scale_CIS_actd'
    'scale_CIS_cisoe'
    'scale_CIS_long'
    'scale_CIS_shp1oe'
    'scale_CIS_socs3oe'
    'scale_SHP1_shp1oe'
    'scale_SOCS3RNA_foldA'
    'scale_SOCS3RNA_foldB'
    'scale_SOCS3RNA_foldC'
    'scale_SOCS3_cisoe'
    'scale_SOCS3_long'
    'scale_SOCS3_socs3oe'
    'scale_pEpoR_actd'
    'scale_pEpoR_cisoe'
    'scale_pEpoR_cisoe_pepor'
    'scale_pEpoR_dr30'
    'scale_pEpoR_dr7'
    'scale_pEpoR_fine'
    'scale_pEpoR_long'
    'scale_pEpoR_shp1oe'
    'scale_pEpoR_socs3oe'
    'scale_pJAK2_actd'
    'scale_pJAK2_cisoe'
    'scale_pJAK2_dr30'
    'scale_pJAK2_dr7'
    'scale_pJAK2_fine'
    'scale_pJAK2_long'
    'scale_pJAK2_shp1oe'
    'scale_pJAK2_socs3oe'
    'scale_pSTAT5_actd'
    'scale_pSTAT5_cisoe'
    'scale_pSTAT5_dr10'
    'scale_pSTAT5_long'
    'scale_pSTAT5_shp1oe'
    'scale_pSTAT5_socs3oe'
    'scale_tSTAT5_actd'
    'scale_tSTAT5_long'
    'scale_tSTAT5_shp1oe'
    'sd_CIS_abs'
    'sd_CIS_au'
    'sd_JAK2EpoR_au'
    'sd_RNA_fold'
    'sd_SHP1_abs'
    'sd_SHP1_au'
    'sd_SOCS3_abs'
    'sd_SOCS3_au'
    'sd_STAT5_abs'
    'sd_STAT5_au'
    'sd_pSTAT5_rel'
    };
parameters.number = nPar; % 112
parameters.min = minPar(1:nPar,1);
parameters.max = maxPar(1:nPar,1);
parameters.guess = par0(1:nPar,1:nStarts);

end
