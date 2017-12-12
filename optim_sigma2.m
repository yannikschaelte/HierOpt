function [varargout] = optim_sigma2(varargin)
% This function computes the optimal the variances of the
% measurements of a model
%
% USAGE:
% [sigma2] = optim_sigma2(iy,simulation,D,distr,options,c)
% [sigma2,dsigma2] = optim_sigma2(iy,simulation,D,distr,options,c,dc)
%
% Parameters
%  iy: index of observable for which the proportionality factor and the variane
%      should be computed
%  for other inputs see nlLH_fgh.m and optim_c.m
%
% Return values:
%   sigma2: optimal variance corresponding to observable i
%   dsigma2: (1 x 1 x n_theta x n_r) matrix with the derivatives of the variances
%   ddsigma2: (n_theta x 1 x n_theta x n_r) matrix that contains the FIM of the
%                                       variances
%% INPUT ASSIGNMENT
iy = varargin{1};
simulation = varargin{2};
D = varargin{3};
distr = varargin{4};
options = varargin{5};
c = varargin{6};
if nargout >= 2
    dc = varargin{7};
end

%% DETERMINE SCENARIO
if(strcmp(options.obs(iy(1)).variance,'multiple') && ...
        strcmp(options.obs(iy(1)).proportionality,'single'))
    warning(['In options.(',num2str(iy(1)),...
        ') you combined variance:multiple and proportionality:single which is not valid.']);
end

%% INITIALIZATION OF DIMENSIONS
n_e = size(D,2); %number of experiments
n_y = size(D(1).my,2); %number of observables
n_r = size(D(1).my,3); %number of replicates
if nargout > 1
    n_theta = size(simulation(1).sy,3); %number of parameters
end

%% COMPUTATION OF SIMGA2, DSIGMA2 AND DDSIGMA2
sigma2 = zeros(1,1,n_r); %variances
if nargout > 1
    dsigma2 = zeros(1,1,n_theta,n_r); %derivatives of the variances
    if nargout > 2
        ddsigma2 = zeros(n_theta,1,n_theta,n_r); %second order derivatives of the variances
    end
end

try
    switch options.obs(iy(1)).variance
        case 'multiple'
            mult_fact_ir = zeros(1,1,n_r);
            sigma2ir = zeros(1,1,n_r);
            if nargout > 1
                dsigma2ir = zeros(1,1,n_theta,n_r);
                if nargout > 2
                    ddsigma2ir = zeros(n_theta,1,n_theta,n_r);
                end
            end
            
            for j = 1:n_e %loop over all experiments
                mult_fact_ir = mult_fact_ir + sum(sum(~isnan(D(j).my(:,iy,:)),1),2);
                switch distr
                    case {'log-laplace','log-normal'}
                        y_ch = bsxfun(@minus,log(D(j).my(:,iy,:)),log(bsxfun(@times,c(:,iy,:,j),simulation(j).y(:,iy))));
                        if nargout > 1
                            dcy_cdy = permute(bsxfun(@ldivide,bsxfun(@times,c(:,iy,:,j),simulation(j).y),...
                                permute(bsxfun(@plus,bsxfun(@times,dc(:,iy,:,:,j),simulation(j).y),...
                                bsxfun(@times,permute(c(:,iy,:,j),[1,2,4,3]),simulation(j).sy)),[1,2,4,3])),[1,2,4,3]);
                        end
                    case {'log10-normal'}
                        y_ch = bsxfun(@minus,log10(D(j).my(:,iy,:)),log10(bsxfun(@times,c(:,iy,:,j),simulation(j).y(:,iy))));
                    case {'laplace','normal'}
                        y_ch = bsxfun(@minus,D(j).my(:,iy,:),bsxfun(@times,c(:,iy,:,j),simulation(j).y(:,iy)));
                        if nargout > 1
                            dcy_cdy = bsxfun(@plus,bsxfun(@times,dc(:,iy,:,:,j),simulation(j).y(:,iy)),...
                                bsxfun(@times,permute(c(:,iy,:,j),[1,2,4,3]),simulation(j).sy(:,iy,:)));
                        end
                end
                
                
                switch distr
                    case {'normal','log-normal'}
                        %computing the optimal variances sigma2_ir
                        sigma2ir = sigma2ir + sum(nansum(bsxfun(@power,y_ch,2),1),2);
                        
                        if nargout > 1
                            %computing the derivatives for the sigma2_ir with respect to theta
                            dsigma2ir = dsigma2ir - sum(nansum(bsxfun(@times,permute(2*y_ch,[1,2,4,3]),dcy_cdy),1),2);
                            
                            if nargout > 2
                                %computing the FIM of the variances with respect to theta
                                error('todo')
                                ddsigma2ir = ddsigma2ir + permute(nansum(...
                                    bsxfun(@times,permute(2*dcy_cdy,[1,2,5,4,3]),dcy_cdy),1),[5,2,3,4,1]);
                            end
                        end
                    case {'laplace','log-laplace'}
                        %computing the optimal sigma_ir
                        sigma2ir = sigma2ir + sum(nansum(abs(y_ch),1),2);
                        if nargout > 1
                            %computing the derivatives for the sigma_ir with respect to theta
                            dsigma2ir = dsigma2ir - sum(nansum(bsxfun(@times,permute(sign(y_ch),[1,2,4,3]),dcy_cdy),1),2);
                        end
                end
            end
            
            sigma2 = bsxfun(@rdivide,sigma2ir,mult_fact_ir);
            if nargout > 1
                dsigma2 = bsxfun(@rdivide,dsigma2ir,permute(mult_fact_ir,[1,2,4,3]));
                if nargout > 2
                    ddsigma2 = bsxfun(@rdivide,ddsigma2ir,permute(mult_fact_ir,[1,2,4,3]));
                end
            end
        case 'single'
            mult_fact_i = zeros(1,1);
            sigma2i = zeros(1,1);
            if nargout > 1
                dsigma2i = zeros(1,1,n_theta);
                if nargout > 2
                    ddsigma2i = zeros(n_theta,1,n_theta);
                end
            end
            
            for j = 1:n_e
                mult_fact_i = mult_fact_i + sum(sum(sum(~isnan(D(j).my(:,iy,:)),1),3));
                
                switch distr
                    case {'log-laplace','log-normal'}
                        y_ch = bsxfun(@minus,log(D(j).my(:,iy,:)),log(bsxfun(@times,c(:,iy,:,j),simulation(j).y(:,iy))));
                        if nargout > 1
                            dcy_cdy = permute(bsxfun(@ldivide,bsxfun(@times,c(:,iy,:,j),simulation(j).y(:,iy)),...
                                permute(bsxfun(@plus,bsxfun(@times,dc(:,iy,:,:,j),simulation(j).y(:,iy)),...
                                bsxfun(@times,permute(c(:,iy,:,j),[1,2,4,3]),simulation(j).sy(:,iy,:))),[1,2,4,3])),[1,2,4,3]);
                        end
                    case {'log10-normal'}
                        y_ch = bsxfun(@minus,log10(D(j).my(:,iy,:)),log10(bsxfun(@times,c(:,iy,:,j),simulation(j).y(:,iy))));
                    case {'laplace','normal'}
                        y_ch = bsxfun(@minus,D(j).my(:,iy,:),bsxfun(@times,c(:,iy,:,j),simulation(j).y(:,iy)));
                        if nargout > 1
                            dcy_cdy = bsxfun(@plus,bsxfun(@times,dc(:,iy,:,:,j),simulation(j).y(:,iy)),...
                                bsxfun(@times,permute(c(:,iy,:,j),[1,2,4,3]),simulation(j).sy(:,iy,:)));
                        end
                end
                
                
                switch distr
                    case {'normal','log-normal','log10-normal'}
                        %computing the optimal variances sigma2_i
                        sigma2i = sigma2i + sum(sum(nansum(bsxfun(@power,y_ch,2),1),3));
                        
                        if nargout > 1
                            %computing the derivatives for the sigma2_i with respect to theta
                            dsigma2i = dsigma2i - sum(sum(nansum(bsxfun(@times,permute(2*y_ch,[1,2,4,3]),dcy_cdy),...
                                1),2),4);
                            
                            if nargout > 2
                                %computing the FIM of the variances with respect to theta
                                ddsigma2i = ddsigma2i + permute(sum(nansum(...
                                    bsxfun(@times,permute(2*dcy_cdy,[1,2,5,4,3]),dcy_cdy),1),4),[5,2,3,4,1]);
                            end
                        end
                    case {'laplace','log-laplace'}
                        %computing the optimal variances sigma_i
                        sigma2i = sigma2i + sum(sum(sum(nansum(abs(y_ch),1),3)));
                        
                        if nargout > 1
                            %computing the derivatives for the sigma_i with respect to theta
                            dsigma2i = dsigma2i - sum(sum(nansum(bsxfun(@times,permute(sign(y_ch),[1,2,4,3]),dcy_cdy),1),2),4);
                        end
                end
            end
            sigma2i = bsxfun(@rdivide,sigma2i,mult_fact_i);
            sigma2 = bsxfun(@times,ones(1,1,n_r),sigma2i);
            
            if nargout > 1
                dsigma2i = bsxfun(@rdivide,dsigma2i,mult_fact_i);
                dsigma2 = bsxfun(@times,ones(1,1,n_theta,n_r),dsigma2i);
                
                if nargout > 2
                    ddsigma2i = bsxfun(@rdivide,ddsigma2i,mult_fact_i);
                    ddsigma2 = bsxfun(@times,ones(n_theta,1,n_theta,n_r),ddsigma2i);
                end
            end
            
    end
catch
    error(['Field options(',num2str(iy(1)),').variance not valid. Valid options: multiple,single.']);
end

switch nargout
    case {0,1}
        varargout{1} = sigma2;
    case 2
        varargout{1} = sigma2;
        varargout{2} = dsigma2;
    case 3
        varargout{1} = sigma2;
        varargout{2} = dsigma2;
        varargout{3} = ddsigma2;
end

% Note: If for all j there are Nans in D(j).my(:,i,r)
% cir and sigma2ir are NaN aswell.

end

