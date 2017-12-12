function [varargout] = nlLH_fgh(varargin)
% This function computes the value of the negative-log-likelihood function,
% its gradient and its Hessian at theta
%
% USAGE:
% [...]                     = nlLH_fgh(simulation,D,'normal')
% [...]                     = nlLH_fgh(simulation,D,'normal',options,1)
% [nllH]                    = nlLH_fgh(...)
% [nllH,gradnlLH]           = nlLH_fgh(...)
% [nllH,gradnlLH,FIMnlLH]   = nlLH_fgh(...)
%
% Parameters:
%   simulation: (1 x #experiments) struct with fields
%       * y: simulation of the output for the different experiments for a theta
%            in which the values nlLH, gradnlLH, FIMnlLH will be computed
%       * sy: simulation of sy, the sensitivities of the output vector, for the different experiments
%            for a theta in which the values nlLH, gradnlLH, FIMnlLH will be computed
%   D: (1 x #experiments) struct containing data with two fields
%       * t: time points
%       * my: number time points x number observables x number replicates
%   distr: indicates distribution of measurement noise,
%         = 'normal'
%         = 'laplace'
%         = 'log-normal'
%         = 'log-laplace'
%   options.exp_groups:
%               * variance: struct with each struct containing the indices
%                   of the experiments that share a variance parameter
%               * proportionality: struct with each struct containing the indices
%                   of the experiments that share a proportionality parameter
%   options.obs_groups:
%               * variance: struct with each struct containing the indices
%                   of the observables that share a variance parameter
%               * proportionality: struct with each struct containing the indices
%                   of the observables that share a proportionality parameter
%    (the defaults for the above options are that all parameters are shared
%    across observables and experiments)
%   options.obs: (1 x #observables) struct with fields
%               * variance:
%                           = 'single': the variance is computed for all replicate together
%                           = 'multiple': the variance is computed for each replicate separately
%               * proportionality:
%                           = 'single': proportionality factor is computed
%                                   for all replicates together
%                           = 'multiple': proportionality factor is
%                                   computed for each replicate separately
%                           = 'absolute': proportionality factor is 1
%   save_pv: indicates whether the proportionality factors and variances are
%            saved in a .mat file (1) or not (0)
%
% Return values:
%   nlLh: value of the negative-log-likelihood function in theta
%   gradnlLH: value of gradient of the negative-log-likelihood function
%                 in theta
%   FIMnlLH: approximation of the Hessian of the negative-log-likelihood function in theta
%            using the FIM
%
% Note: all oberservables/experiments that share a proportionality parameter
%       need to also share the variance parameter!
%
%% INITIALIZATION
simulation = varargin{1};
D = varargin{2};
distr = varargin{3};
options = varargin{4};
if nargin >= 5
    save_pv = varargin{5};
else
    save_pv = 0;
end

n_e = size(D,2); %number of experiments
n_y = size(D(1).my,2); %number of observables
n_r = size(D(1).my,3); %number of replicates

if ~isfield(options,'exp_groups')
    options.exp_groups.proportionality{1} = [1:n_e];
    options.exp_groups.variance{1} = [1:n_e];
else
    if ~isfield(options.exp_groups,'proportionality')
        options.exp_groups.proportionality{1} = [1:n_e];
    end
    if ~isfield(options.exp_groups,'variance')
        options.exp_groups.variance{1} = [1:n_e];
    end
end

if ~isfield(options,'obs_groups')
    options.obs_groups.proportionality{1} = [1:n_y];
    options.obs_groups.variance{1} = [1:n_y];
else
    if ~isfield(options.obs_groups,'proportionality')
        options.obs_groups.proportionality{1} = [1:n_y];
    end
    if ~isfield(options.obs_groups,'variance')
        options.obs_groups.variance{1} = [1:n_y];
    end
end

if ~isfield(options,'lsqnonlin')
    options.lsqnonlin = false;
end

if nargout > 1
    n_theta = size(simulation(1).sy,3);%number of parameters
end

nlLH = 0;
if nargout > 1
    gradnlLH = zeros(n_theta,1);
    if nargout >  2
        if(strcmp(distr,'laplace'))
            error('No user supplied Fisher information matrix available for the Laplace distribution.');
        else
            FIMnlLH = zeros(n_theta,n_theta);
        end
    end
end

% Initialization for for proportionality and variance parameters
c = zeros(1,n_y,n_r,n_e); % vector including proportionality factors
sigma2 = zeros(1,n_y,n_r,n_e); % vector including variances

if nargout > 1
    dc = zeros(1,n_y,n_theta,n_r,n_e);
end

% if ~checkValidityOptions(D,options)
%     error(['options not valid, all oberservables/experiments that share a' ...
%         'proportionality parameter need to also share the variance parameter'])
% end

%% OPTIMAL VALUES FOR THE PROPORTIONALITY FACTORS
for ie = 1:numel(options.exp_groups.proportionality)
    for iy = 1:numel(options.obs_groups.proportionality)
        ind_y = options.obs_groups.proportionality{iy};
        ind_e = options.exp_groups.proportionality{ie};

        switch distr
            case {'normal','log-normal','log10-normal'}
                tempc = optim_c(ind_y,simulation(ind_e),D(ind_e),distr,options);
                for ir = 1:n_r
                    c(:,ind_y,ir,ind_e) = tempc(:,:,ir);
                end
            case {'laplace','log-laplace'}
                if nargout > 1
                    [tempc,tempdc] = ...
                        optim_c(ind_y,simulation(ind_e),D(ind_e),distr,options);
                    for ir = 1:n_r
                        c(:,ind_y,ir,ind_e) = tempc(:,:,ir);
                        for iiy = 1:numel(ind_y)
                            for iie = 1:numel(ind_e)
                                dc(:,ind_y(iiy),:,ir,ind_e(iie)) = tempdc(:,:,:,ir);
                            end
                        end
                    end
                else
                    tempc = optim_c(ind_y,simulation(ind_e),D(ind_e),distr,options);
                    for ir = 1:n_r
                        c(:,ind_y,ir,ind_e) = tempc(:,:,ir);
                    end
                end
        end
    end
end

%% OPTIMAL VALUES FOR THE VARIANCES
for ie = 1:numel(options.exp_groups.variance)
    for iy = 1:numel(options.obs_groups.variance)
        ind_y = options.obs_groups.variance{iy};
        ind_e = options.exp_groups.variance{ie};
        [tempsigma2] = ...
            optim_sigma2(ind_y,simulation(ind_e),D(ind_e),distr,options,c(:,:,:,ind_e));
        for ir = 1:n_r
            sigma2(:,ind_y,ir,ind_e)=tempsigma2(:,:,ir);
        end
        
    end
end

%% SAVING PROPORTIONALITY FACTORS AND VARIANCES
if(save_pv == 1)
    save('analytical_results.mat','c','sigma2');
end

%% COMPUTING nlLH, gradnlLH, FIMnlLH
if options.lsqnonlin
    res = [];
    sres = [];
end
switch distr
    case {'normal','log-normal','log10-normal'}
        for j = 1:n_e
            sigma2_j = sigma2(:,:,:,j);
            c_j = c(:,:,:,j);
            switch distr
                case 'normal'
                    y_ch = bsxfun(@minus,D(j).my,bsxfun(@times,c_j,simulation(j).y));
                    if nargout > 1
                        dy_ch = -bsxfun(@times,permute(c_j,[1,2,4,3]),repmat(simulation(j).sy,[1,1,1,n_r]));
                    end
                case 'log-normal'
                    y_ch = bsxfun(@minus,log(D(j).my),log(bsxfun(@times,c_j,simulation(j).y)));
                    if nargout > 1
                        dy_ch = -bsxfun(@ldivide,simulation(j).y,repmat(simulation(j).sy,[1,1,1,n_r]));
                    end
                    if options.lsqnonlin
                        error('check this options, is dsigma2 needed?')
                        for ir = 1:size(D(j).my,3)
                            temp = bsxfun(@rdivide,y_ch(:,:,ir),sqrt(sigma2(:,:,ir)));
                            res = [res;temp(:)];
                            temp_res_err = sqrt(bsxfun(@times,~isnan(D(j).my(:,:,ir)),log(2*pi*sigma2(:,:,ir))+50));
                            res_err = temp_res_err(:);
                            ind = find(res_err>0);
                            res = [res;res_err(ind)];
                            if nargout > 1
                                stemp = bsxfun(@minus,-bsxfun(@times,temp./sqrt(sigma2(:,:,ir)),...
                                    bsxfun(@rdivide,dsigma2(:,:,:,ir),2*sqrt(sigma2(:,:,ir)))),...
                                    bsxfun(@times,1./(bsxfun(@times,sqrt(sigma2(:,:,ir)),...
                                    bsxfun(@times,c_j,simulation(j).y))),dcy_cdy));
                                stemp = reshape(stemp,numel(temp),size(simulation(j).sy,3));
                                sres = [sres; stemp];
                                temp_sres_err = bsxfun(@times,bsxfun(@rdivide,1./temp_res_err,2*sigma2(:,:,ir)),...
                                    dsigma2(:,:,:,ir));
                                temp_sres_err = reshape( temp_sres_err,numel(temp_res_err),size(simulation(j).sy,3));
                                sres = [sres; temp_sres_err(ind,:)];
                            end
                        end
                    end
                case 'log10-normal'
                    y_ch = bsxfun(@minus,log10(D(j).my),log10(bsxfun(@times,c_j,simulation(j).y)));
                    if nargout > 1
                        dy_ch = -bsxfun(@ldivide,simulation(j).y*log(10),repmat(simulation(j).sy,[1,1,1,n_r]));
                    end
            end
            if ~options.lsqnonlin
                nlLH = nlLH + sum(sum(nansum(bsxfun(@times,~isnan(D(j).my),log(2*pi*sigma2_j))+...
                    bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2_j),1),3),2);
                
                if nargout > 1
                    %value of the gradient of the neg-log-likelihood function
                    gradnlLH = gradnlLH + 2*permute(sum(sum(nansum(...
                        bsxfun(@times,permute(...
                        bsxfun(@rdivide,y_ch,sigma2_j),[1,2,4,3]),...
                        dy_ch),1),4),2),[3,2,1]);
                    
                    if nargout > 2
                        %Fisher information matrix
                        error('to do')
                        ddsigma2_j = ddsigma2(:,:,:,:,:,j);
                        FIMnlLH = FIMnlLH +...
                            permute(sum(sum(nansum(bsxfun(@rdivide,...
                            bsxfun(@times,permute(~isnan(D(j).my),[1,2,4,3]),-bsxfun(@rdivide,bsxfun(@times,...
                            permute(dsigma2_j,[1,2,5,4,3]),dsigma2_j),permute(sigma2_j,[1,2,4,3]))+...
                            permute(ddsigma2_j,[5,2,3,4,1]))+...
                            bsxfun(@times,permute(bsxfun(@times,2*c_j,simulation(j).y),[1,2,4,3]),...
                            bsxfun(@times,permute(dc_j,[1,2,5,4,3]),simulation(j).sy)+...
                            bsxfun(@times,dc_j,permute(simulation(j).sy,[1,2,5,4,3])))+...
                            bsxfun(@times,2*bsxfun(@power,simulation(j).y,2),...
                            bsxfun(@times,permute(dc_j,[1,2,5,4,3]),dc_j))+...
                            bsxfun(@times,permute(2*bsxfun(@power,c_j,2),[1,2,4,3]),...
                            bsxfun(@times,permute(simulation(j).sy,[1,2,5,4,3]),simulation(j).sy)),...
                            sigma2),1),4),2),[5,3,2,4,1]);
                    end
                end
            end
        end
        nlLH = 0.5*nlLH;
        if nargout > 1
            gradnlLH = 0.5*gradnlLH;
            if nargout > 2
                FIMnlLH = 0.5*FIMnlLH;
            end
        end
    case {'laplace','log-laplace'}
        for j = 1:n_e
            sigma2_j = sigma2(:,:,:,j);
            c_j = c(:,:,:,j);
            if nargout > 1
                dc_j = dc(:,:,:,:,j);
            end
            if strcmp(distr,'laplace')
                y_ch = bsxfun(@minus,D(j).my,bsxfun(@times,c_j,simulation(j).y));
                if nargout > 1
                    dcy_cdy = bsxfun(@times,dc_j,simulation(j).y) +...
                        bsxfun(@times,permute(c_j,[1,2,4,3]),simulation(j).sy);
                end
            else
                y_ch = bsxfun(@minus,log(D(j).my),log(bsxfun(@times,c_j,simulation(j).y)));
                if nargout > 1
                    dcy_cdy = permute(bsxfun(@ldivide,bsxfun(@times,c_j,simulation(j).y),...
                        permute(bsxfun(@plus,bsxfun(@times,dc_j,simulation(j).y),...
                        bsxfun(@times,permute(c_j,[1,2,4,3]),simulation(j).sy)),[1,2,4,3])),[1,2,4,3]);
                end
            end
            nlLH = nlLH + sum(sum(nansum(bsxfun(@times,~isnan(D(j).my),log(2*sigma2_j))+...
                bsxfun(@rdivide,abs(y_ch),sigma2_j),1),3),2);
            if nargout > 1
                gradnlLH = gradnlLH - permute(sum(sum(nansum(...
                    bsxfun(@times,permute(bsxfun(@rdivide,...
                    sign(y_ch),sigma2_j),[1,2,4,3]),dcy_cdy)...
                    ,1),4),2),[3,2,1]);
            end
        end
end
if isfield(options,'lsqnonlin') && options.lsqnonlin
    switch nargout
        case{0,1}
            varargout{1} = res;
        case 2
            varargout{1} = res;
            varargout{2} = sres;
    end
else
    switch nargout
        case{0,1}
            varargout{1} = nlLH;
        case 2
            varargout{1} = nlLH;
            varargout{2} = gradnlLH;
        case 3
            varargout{1} = nlLH;
            varargout{2} = gradnlLH;
            varargout{3} = FIMnlLH;
    end
end
% varargout{1} = simulation(1).y(1,1);
% if nargout >=2
% varargout{2} = squeeze(simulation(1).sy(1,1,:));
% end
end

