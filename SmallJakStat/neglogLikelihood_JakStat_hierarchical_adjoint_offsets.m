function [varargout] = neglogLikelihood_JakStat_hierarchical_adjoint(varargin)
% compute the loglikelihood, gradient and Hessian for RafMekErk data

%INPUT:
%xi: parameter vector
%kappa: constants vector
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
    theta = varargin{1};
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


    n_e = size(D,2); %number of experiments
    n_y = size(D(1).my,2); %number of observables
    n_r = size(D(1).my,3); %number of replicates
    
    amioptions = amioption();
    amioptions.sensi_meth = 'adjoint'; % 2
    
    %% negative log-likelihood
    % first run: simulate absolute observables
    amioptions.sensi = 0; % do not compute derivatives
    
    sim = struct([]);
    for j = 1:n_e
        % in every experiment, the timepoints and conditions are the same
        % over all observables and replicates
        tout_j = D(j).t;
        kappa_j = [ kappa(j,:) 0 0 1 1 ]; % default values for scalings
        data_j = []; % no data
        sol = simulate_JakStat_hierarchical_adjoint_offsets(tout_j,theta,kappa_j,data_j,amioptions);
        if sol.status < 0
            error(['failed to integrate ODE for experiment ' num2str(j)])
        end
        sim(j).y = sol.y;
    end
    
    % scalings
    options.obs(1).variance = 'multiple';
    options.obs(1).proportionality = 'multiple';
    options.obs(2).variance = 'multiple';
    options.obs(2).proportionality = 'multiple';
    options.obs(3).variance = 'multiple';
    options.obs(3).proportionality = 'absolute';
    options.obs_groups.variance = {1,2,3};
    options.obs_groups.proportionality = {1,2,3};
    
    [c,sigma2,c_by_y,sigma2_by_y] = scalings(sim,D,options);
    
    % nllh
    nlLH = 0;
    for ie = 1:n_e
        sigma2_e = sigma2(:,:,:,ie);
        c_e = c(:,:,:,ie);
        b_e = zeros(size(c_e));
        y_ch = bsxfun(@minus,D(ie).my,bsxfun(@times,c_e,sim(ie).y)+b_e);
        nlLH = nlLH + sum(sum(nansum(bsxfun(@times,~isnan(D(ie).my),log(2*pi*sigma2_e))+...
                    bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2_e),1),3),2);
    end
    nlLH = 0.5*nlLH;
    varargout{1} = nlLH;
    
    %% derivatives
    % second run: simulate with correct scalings and adjoints
    % if we want to compute derivatives, run the simulation again with the
    % correct scalings and compute the derivative of J using adjoints
    if nargout > 1
        amioptions.sensi = 1;
        amioptions.sensi_meth = 'adjoint';
        snlLH = 0;
        nlLH2 = 0;
        for ie = 1:n_e
            for ir = 1:n_r
                tout_e = D(ie).t;
                
                c_re = c_by_y(:,:,ir,ie);
                c_re = reshape(c_re,1,[]);
                
%                 b_re = b_by_y(:,:,ir,ie);
                b_re = zeros(size(c_re));
                b_re = reshape(b_re,1,[]);
                
                sigma2_re = sigma2_by_y(:,:,ir,ie);
                sigma2_re = reshape(sigma2_re,1,[]);
     
                sigma2_re = repmat(sigma2_re,length(tout_e),1);
                
                kappa_re = [ kappa(ie,:) b_re c_re ]; % default values for scalings
                
                my_re = D(ie).my(:,:,ir); % data
                data_re.t = tout_e;
                data_re.Y = my_re;
                data_re.Sigma_Y = sqrt(sigma2_re);
                data_re.condition = kappa_re;
                data_re = amidata(data_re);
                
                sol = simulate_JakStat_hierarchical_adjoint_offsets(tout_e,theta,kappa_re,data_re,amioptions);
                if sol.status < 0
                    error(['failed to integrate ODE for experiment ' num2str(ie)])
                end
                nlLH2 = nlLH2 - sol.llh;
                snlLH = snlLH - sol.sllh;
            end
        end
%         fprintf('%.15f %.15f\n',nlLH,nlLH2);
        varargout{2} = snlLH;
    end
