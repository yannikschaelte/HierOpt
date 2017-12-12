function [varargout] = optim_c(iy,simulation,D,distr,options)
% This function computes the optimal proportionality factors of the
% measurements of a model
%
% USAGE:
% [...]  = optim_c(iy,simulation,D,distr,options)
% [c,dc] = optim_c(...)
%
% Parameters:
%  iy: observable index for which the proportionality factor and the
%      variance should be computed
%  for the other inputs see nlLH_fgh.m
%
% Return values:
%   c: the optimal proportionality factor for observable i
%   dc: (1 x 1 x n_theta x n_r) matrix with the derivatives of the proportionality factors

%% CHECK SCENARIO
for iiy = 1:numel(iy)
    if(strcmp(options.obs(iy(iiy)).variance,'multiple') && ...
            strcmp(options.obs(iy(iiy)).proportionality,'single'))
        error(['In options.(',num2str(iy(iiy)),...
            ') you combined variance:multiple and proportionality:single which is not valid.']);
    end
end
% all data points that are used for the calculation of a proportionality
% factors need to have the same variance parameter!

proportionality = options.obs(iy(1)).proportionality;
variance = options.obs(iy(1)).variance;
for iiy = 1:numel(iy)
    if ~strcmp(options.obs(iy(iiy)).variance,variance) || ...
            ~strcmp(options.obs(iy(iiy)).proportionality,proportionality)
        if ~strcmp(options.obs(iy(iiy)).proportionality,'absolute')
            error('different assumptions')
        end
    end
end

%% INITIALIZATION OF DIMENSIONS
n_e = size(D,2); %number of experiments
n_r = size(D(1).my,3); %number of replicates
if nargout > 1
    n_theta = size(simulation(1).sy,3); %number of parameters
end
%% COMPUTATION OF C AND DC
c = zeros(1,1,n_r); %proportionality factors
if nargout > 1
    dc = zeros(1,1,n_theta,n_r); %derivatives of the proportionality factors
end

try
    switch proportionality
        case 'absolute'
            c = ones(1,1,n_r);
            if nargout > 1
                dc = zeros(1,1,n_theta,n_r);
            end
        case 'multiple'
            switch distr
                case {'normal'}
                    cir_z = zeros(1,1,n_r);
                    cir_n = zeros(1,1,n_r);
                    if nargout > 1
                        dcir_naz = zeros(1,1,n_theta,n_r);
                        dcir_zan = zeros(1,1,n_theta,n_r);
                    end
                    
                    for j = 1:n_e %loop over all experiments
                        %calculating the optimal proportionality factors c_ir
                        cir_z = cir_z + nansum(nansum(bsxfun(@times,D(j).my(:,iy,:),simulation(j).y(:,iy)),1),2);
                        cir_n = cir_n + sum(sum(bsxfun(@power,bsxfun(@times,~isnan(D(j).my(:,iy,:)),...
                            simulation(j).y(:,iy)),2),1),2);
                        
                        if nargout > 1
                            %computing the derivatives of c_ir with respect to theta
                            dcir_naz = dcir_naz + sum(nansum(bsxfun(@times,simulation(j).sy(:,iy,:),...
                                permute(D(j).my(:,iy,:),[1,2,4,3])),1),2);
                            dcir_zan = dcir_zan + sum(sum(bsxfun(@times,...
                                permute(bsxfun(@times,~isnan(D(j).my(:,iy,:)),...
                                2*simulation(j).y(:,iy)),[1,2,4,3]),simulation(j).sy(:,iy,:)),1),2);
                        end
                    end
                    c = bsxfun(@rdivide,cir_z,cir_n);
                    if nargout > 1
                        dc = bsxfun(@minus,bsxfun(@rdivide,dcir_naz,permute(cir_n,[1,2,4,3])),...
                            bsxfun(@rdivide,bsxfun(@times,permute(cir_z,[1,2,4,3]),dcir_zan),...
                            permute(bsxfun(@power,cir_n,2),[1,2,4,3])));
                    end
                case {'log-normal'}
                    error('not yet implemented')
                    
                case {'laplace'}
                    for ir = 1:n_r
                        candidates = [];
                        grad_candidates = [];
                        for j = 1:n_e
                            for iiy = 1:numel(iy)
                                for it = 1:size(D(j).my,1)
                                    if ~isnan(D(j).my(it,iy(iiy),ir))
                                        candidates = [candidates;bsxfun(@rdivide,...
                                            D(j).my(it,iy(iiy),ir),simulation(j).y(it,iy(iiy)))];
                                        if nargout > 1
                                            grad_candidates = [grad_candidates;-bsxfun(@rdivide,...
                                                simulation(j).sy(it,iy(iiy),:),simulation(j).y(it,iy(iiy)))];
                                        end
                                    end
                                end
                            end
                        end
                        [candidates,I] = sort(candidates); %candidates for c_ir, should be n_e*n_y*n_t
                        middle = (candidates(1:end-1,:,:)+candidates(2:end,:,:))/2;
                        dJdc = zeros(size(middle));
                        
                        for cand = 1:size(dJdc,1) %calculating for all candidates (sigma*) derivative of J with respect to c_ir
                            for j = 1:n_e
                                for iiy = 1:numel(iy)
                                    for it = 1:size(D(j).my,1)
                                        if ~isnan(D(j).my(it,iy(iiy),ir))
                                            yh_c = bsxfun(@minus,bsxfun(@rdivide,...
                                                D(j).my(it,iy(iiy),ir),simulation(j).y(it,iy(iiy))),middle(cand,:,:));
                                            dJdc(cand,1) = dJdc(cand,1) - nansum(bsxfun(@times,...
                                                abs(simulation(j).y(it,iy(iiy))),sign(yh_c)),1);
                                            
                                        end
                                    end
                                end
                            end
                        end
                        c_opt_r = find(dJdc==0);
                        if isempty(c_opt_r)
                            c_opt_r = find(bsxfun(@times,dJdc(1:end-1,1)<= 0,dJdc(2:end,1,1)>=0))+1;
                        end
                        if ~isempty(c_opt_r)
                            c(1,1,ir) = candidates(c_opt_r);
                        else
                            if dJdc(1,1) > 0
                                c_opt_r = 1;
                            else
                                c_opt_r = size(dJdc,1); %never happens
                            end
                            c(1,1,ir) = candidates(c_opt_r);
                        end
                        
                        if nargout > 1
                            grad_candidates = grad_candidates(I,:);
                            dc_i = c(1,1,ir)*squeeze(grad_candidates(c_opt_r,:,:));
                            dc(1,1,:,ir) = dc_i';
                        end
                    end
            end
        case 'single'
            switch distr
                case {'log-normal','log10-normal'}
                    logmy = zeros(1,1);
                    logy = zeros(1,1);
                    if nargout > 1
                        dci_part = zeros(1,1,n_theta);
                    end
                    multfact = 0;
                    for j = 1:n_e
                        %calculating the optimal proportionality factors c_ir
                        logmy = logmy + sum(sum(sum(nansum(log(D(j).my(:,iy,:))))));
                        logy  = logy + sum(sum(sum((bsxfun(@times,~isnan(D(j).my(:,iy,:)),...
                            log(simulation(j).y(:,iy)))))));
                        multfact = multfact + sum(sum(sum(~isnan(D(j).my(:,iy,:)))));
                        if nargout > 1
                            %computing the derivatives of c_ir with respect to theta
                            tempD = double(~isnan(D(j).my(:,iy,:)));
                            tempD(tempD==0) = NaN;
                            dci_part = dci_part + sum(sum(nansum(bsxfun(@ldivide,permute(bsxfun(@times,...
                                tempD,simulation(j).y(:,iy)),[1,2,4,3]),...
                                simulation(j).sy(:,iy,:)),1),2),4);
                        end
                    end
                    c = bsxfun(@times,ones(1,1,n_r),exp((logmy-logy)/multfact));
                    if nargout > 1
                        dc = bsxfun(@times,ones(1,1,n_theta,n_r),bsxfun(@times,c(1)/multfact,... %check whether valid
                            -dci_part));
                    end
                case {'normal'}
                    ci_z = zeros(1,1);
                    ci_n = zeros(1,1);
                    if nargout > 1
                        dci_naz = zeros(1,1,n_theta);
                        dci_zan = zeros(1,1,n_theta);
                    end
                    
                    for j = 1:n_e
                        %calculating the optimal proportionality factors c_ir
                        ci_z = ci_z + sum(sum(nansum(bsxfun(@times,D(j).my(:,iy,:),simulation(j).y(:,iy)),1),3));
                        ci_n = ci_n + sum(sum(sum(bsxfun(@power,bsxfun(@times,~isnan(D(j).my(:,iy,:)),...
                            simulation(j).y(:,iy)),2),1),3));
                        
                        if nargout > 1
                            %computing the derivatives of c_ir with respect to theta
                            dci_naz = dci_naz + sum(sum(nansum(bsxfun(@times,permute(D(j).my(:,iy,:),[1,2,4,3]),...
                                simulation(j).sy(:,iy,:)),1),2),4);
                            dci_zan = dci_zan + sum(sum(sum(bsxfun(@times,permute(bsxfun(@times,...
                                ~isnan(D(j).my(:,iy,:)),2*simulation(j).y(:,iy)),[1,2,4,3]),...
                                simulation(j).sy(:,iy,:)),1),2),4);
                        end
                    end
                    
                    c = bsxfun(@times,ones(1,1,n_r),bsxfun(@rdivide,ci_z,ci_n));
                    if nargout > 1
                        dc = bsxfun(@times,ones(1,1,n_theta,n_r),bsxfun(@minus,bsxfun(@rdivide,dci_naz,ci_n),...
                            bsxfun(@rdivide,bsxfun(@times,ci_z,dci_zan),bsxfun(@power,ci_n,2))));
                    end
                case {'laplace','log-laplace'}
                    candidates = [];
                    grad_candidates = [];
                    for j = 1:n_e
                        for iiy = 1:numel(iy)
                            for it = 1:size(D(j).my,1)
                                for ir = 1:n_r
                                    if ~isnan(D(j).my(it,iy(iiy),ir))
                                        candidates = [candidates;reshape(bsxfun(@rdivide,...
                                            D(j).my(it,iy(iiy),ir),simulation(j).y(it,iy(iiy))),[],1)];
                                        if nargout > 1
                                            grad_candidates = [grad_candidates;repmat(-bsxfun(@rdivide,...
                                                simulation(j).sy(it,iy(iiy),:),simulation(j).y(it,iy(iiy))),[n_r,1])];
                                        end
                                    end
                                end
                            end
                            
                        end
                    end
                    [candidates,I] = sort(candidates); %candidates for c_ir, should be n_e*n_y*n_t
                    middle = (candidates(1:end-1,:,:)+candidates(2:end,:,:))/2;
                    dJdc = zeros(size(middle));
                    
                    for cand = 1:size(dJdc,1)
                        for j = 1:n_e
                            for iiy = 1:numel(iy)
                                for it = 1:size(D(j).my,1)
                                    for ir = 1:n_r
                                        if ~isnan(D(j).my(it,iy(iiy),ir))
                                            yh_c = bsxfun(@minus,bsxfun(@rdivide,...
                                                D(j).my(it,iy(iiy),ir),simulation(j).y(it,iy(iiy))),middle(cand,:,:));
                                            dJdc(cand,1) = dJdc(cand,1) - sum(nansum(bsxfun(@times,...
                                                abs(simulation(j).y(it,iy(iiy))),sign(yh_c)),1),3);
                                            
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if ~isempty(dJdc)
                        c_opt = find(dJdc==0);
                        if isempty(c_opt)
                            c_opt = find(bsxfun(@times,dJdc(1:end-1,1)<= 0,dJdc(2:end,1)>=0))+1;
                        end
                        if isempty(c_opt) && dJdc(1) > 0
                            c_opt = 1;
                        end
                        c = bsxfun(@times,ones(1,1,n_r),candidates(c_opt));
                        
                        if nargout > 1
                            grad_candidates = grad_candidates(I,:);
                            %dc_i = squeeze(grad_candidates(c_opt,:,:));
                            dc_i = c(1,1,1)*squeeze(grad_candidates(c_opt,:,:));
                            for r = 1:n_r
                                dc(1,1,:,r) = dc_i';
                            end
                            %dc = bsxfun(@times,ones(1,1,n_theta,n_r),dc_i');
                        end
                     end
            end
    end
catch
    error(['error for calculation of proportionality of observable ' num2str(iy)]);
end

if(any(c == 0))
    warning(['At least one of the computed proportionality factors for observable is 0.']);
end

switch nargout
    case{0,1}
        varargout{1} = c;
    case 2
        varargout{1} = c;
        varargout{2} = dc;
end

% Note: If for all j there are Nans in D(j).my(:,i,r)
% cir and sigma2ir are NaN as well.

end

