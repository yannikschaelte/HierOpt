function [varargout] = neglogLikelihood_RafMekErk_hierarchical(varargin)
% compute the loglikelihood, gradient and Hessian for RafMekErk data

%INPUT:
%xi: parameter vector
%kappa: konstants vector
%D: data struct of the form D.t D.my
%distr: indicates distribution of measurement noise,
%           * normal
%           * laplace
%save_pv: indicates weather the proportionality factors and variances are
%         saved in a .mat file (1) or not (0)


%OUTPUT
%negative-log-Likelihood
%Gradient of the negative-log-Likelihood

% CHECK/ASSIGN INPUTS:
if nargin >= 4
    xi = varargin{1};
    kappa = varargin{2};
    D = varargin{3};
    distr = varargin{4};
    save_pv = varargin{5};
else
    error('Not enough input arguments!');
end

if nargin == 6
    prior2 = varargin{6};
end

try
    %global errorCount
   n_e = size(D,2);
   if nargout>1
        options_simu.sensi = 1;
    else
        options_simu.sensi = 0;
   end
    %% SIMULATION
    simulation = struct([]);
    for j = 1:n_e %simulations for the different input values for Sora and UO126
            sol = simulate_RafMekErk_hierarchical(D(j).t,xi,kappa(j,:)',[],options_simu);        
        if sol.status < 0 
            error(['failed to integrate ODE for experiment ' num2str(j)])
        end
        
       simulation(j).y = sol.y;
       if nargout > 1
            simulation(j).sy = sol.sy;
       end
    end
    
    %% NEGATIVE-LOG-LIKELIHOOD, GRADIENT, FIM
    
    rep_mode = 'multiple';
    
    %y(1)
    options.obs(1).variance = rep_mode;
    options.obs(1).proportionality = rep_mode;
    %y(2)
    options.obs(2).variance = rep_mode;
    options.obs(2).proportionality = rep_mode;

    options.obs_groups.proportionality{1} = 1;
    options.obs_groups.proportionality{2} = 2;
    options.obs_groups.variance{1} = 1;
    options.obs_groups.variance{2} = 2;


    if nargout > 1
        [nlLH, gradnlLH] = nlLH_fgh(simulation,D,distr,options,save_pv);
        %posterior = likelihood * prior => -log(posterior) = -log(likelihood) - log(prior)
        if(nargin == 6)
            nlLH = nlLH + 0.5*sum(xi.^2)/prior2;
            gradnlLH = gradnlLH + xi/prior2;
        end
    else
        nlLH = nlLH_fgh(simulation,D,distr,options,save_pv);
        %posterior = likelihood * prior => -log(posterior) = -log(likelihood) - log(prior)
        if(nargin == 6)
            nlLH = nlLH + 0.5*sum(xi.^2)/prior2;
        end
    end
    
catch error_thrown
    warning(['Evaluation of likelihood failed. ',error_thrown.message]);
    nlLH = inf; 
    gradnlLH = zeros(length(xi),1);
end

switch nargout
    case{0,1}
        varargout{1} = nlLH;
    case 2
        varargout{1} = nlLH;
        varargout{2} = gradnlLH;
%         varargout{3} = dummy1;
%         varargout{4} = dummy2;
end

end