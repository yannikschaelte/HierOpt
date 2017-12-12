function [parameters,options] = getParameterOptions_RafMekErk(approach)

options.MS = PestoOptions();
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
    'display','iter',...
    'GradObj','on',...
    'MaxIter',5000,...
    'TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxFunEvals',40000,...
    'PrecondBandWidth', inf);
options.MS.n_starts = 10;%500
options.MS.mode = 'text';
options.MS.obj_type = 'negative log-posterior';
options.MS.save = true;

load parameter_guesses_RafMekErk par0

if contains(approach,'hierarchical')
    parameters.name = {'log_{10}(kdf_Raf)','log_{10}(kp_Raf)','log_{10}(kdp_pMek)',...
        'log_{10}(kp_pRaf_Mek)','log_{10}(kdp_pErk)','log_{10}(kp_pMek_Erk)',...
        'log_{10}(K_pErk_inh)','log_{10}(sust_Ras_0)','log_{10}(ts_sust_Ras)',...
        'log_{10}(ts_trans_Ras)','log_{10}(K_Sora)','log_{10}(K_UO)'};
    parameters.guess = par0(1:length(parameters.name),:);
else
    parameters.name = {'log_{10}(kdf_Raf)','log_{10}(kp_Raf)','log_{10}(kdp_pMek)',...
        'log_{10}(kp_pRaf_Mek)','log_{10}(kdp_pErk)','log_{10}(kp_pMek_Erk)',...
        'log_{10}(K_pErk_inh)','log_{10}(sust_Ras_0)','log_{10}(ts_sust_Ras)',...
        'log_{10}(ts_trans_Ras)','log_{10}(K_Sora)','log_{10}(K_UO)',...
        'log_{10}(scale_pMek_20140430_gel1)','log_{10}(scale_pErk_20140430_gel1)',...
        'log_{10}(scale_pMek_20140430_gel2)','log_{10}(scale_pErk_20140430_gel2)',...
        'log_{10}(scale_pMek_20140505_gel1)','log_{10}(scale_pErk_20140505_gel1)',...
        'log_{10}(scale_pMek_20140505_gel2)','log_{10}(scale_pErk_20140505_gel2)',...
        'log_{10}(sigma_pMek_20140430_gel1)','log_{10}(sigma_pErk_20140430_gel1)',...
        'log_{10}(sigma_pMek_20140430_gel2)','log_{10}(sigma_pErk_20140430_gel2)',...
        'log_{10}(sigma_pMek_20140505_gel1)','log_{10}(sigma_pErk_20140505_gel1)',...
        'log_{10}(sigma_pMek_20140505_gel2)','log_{10}(sigma_pErk_20140505_gel2)'
        };
    parameters.guess = par0;
end
parameters.number = length(parameters.name);
parameters.min = -7*ones(parameters.number,1);
parameters.max = 5*ones(parameters.number,1);




