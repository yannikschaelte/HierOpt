function varargout = neglogLikelihood_Bachmann(xi,D,options)
nderiv = nargout-1;
if(nderiv>=1)
    options.ami.sensi = 1;
else
    options.ami.sensi = 0;
end
if ~isfield(options.llh,'lsqnonlin')
    options.llh.lsqnonlin = 0;
end
if ~isfield(options.llh,'reduced_woinit')
    options.llh.reduced_woinit = false;
end
if ~isfield(options.llh,'save_analytical')
    options.llh.save_analytical = false;
end
if options.llh.reduced_woinit
    xi_temp = nan(numel(xi)+2,1);
    xi_temp([1:24,27:end]) = xi;
    xi_temp(25) = 0.59947;
    xi_temp(26) = 1.4269;
    xi = xi_temp;
end
% Simulation of conditions
try
    for cond = 1:numel(D)
        options.ami.x0 = D(cond).init(xi,D(cond).u);
        if options.ami.sensi
            options.ami.sx0 = D(cond).sinit(xi,D(cond).u);
        end
        if options.llh.reduced_woinit
            sol(cond) = simulate_Bachmann_JAKSTAT_red(D(cond).t,xi(1:27),D(cond).u,[],options.ami);
        else
            sol(cond) = simulate_Bachmann_JAKSTAT_D2D(D(cond).t,xi(1:27),D(cond).u,[],options.ami);
        end
        temp_status(cond) = sol(cond).status;
    end
catch
    varargout{1} = NaN;
    if nderiv>=1
        varargout{2} = nan(numel(xi),1);
    end
    warning('simulation failed')
    return;
end
if any(temp_status<0)
    if options.llh.lsqnonlin
        varargout{1} = nan(1082 ,1);
        if nderiv>=1
            varargout{2} = nan(1082 ,numel(xi));
        end
    else
        varargout{1} = nan;
        if nderiv>=1
            varargout{2} = nan(numel(xi),1);
        end
        
    end
    warning('simulation failed');
    return;
end

res = [];
sres = [];
switch options.llh.approach
    case 'hierarchical'
        sol = getSimulation_Bachmann_JAKSTAT_offsetscaling(xi,sol,D,options);
        if options.llh.lsqnonlin
            if nderiv == 0
                res = nlLH_fgh(sol,D,options.llh.distribution,options.llh,options.llh.save_analytical);
            else
                [res,sres] = nlLH_fgh(sol,D,options.llh.distribution,options.llh,options.llh.save_analytical);
            end
        else
            if nderiv == 0
                nlogL = nlLH_fgh(sol,D,options.llh.distribution,options.llh,options.llh.save_analytical);
            else
                [nlogL,dnlogL] = nlLH_fgh(sol,D,options.llh.distribution,options.llh,options.llh.save_analytical);
            end
        end
    case 'standard'
        nlogL = 0;
        dnlogL=zeros(numel(xi),1);
        sol = getSimulation_Bachmann_JAKSTAT_offsetscaling(xi,sol,D,options);
        for cond = 1:numel(D)
            % Map variance parameters
            sigma2 = zeros(1,20,size(D(cond).my,3));
            for r = 1:size(D(cond).my,3)
                if options.llh.original
                    sigma2(1,:,r) = (10.^(2*xi(D(cond).std)));
                else
                    sigma2(1,:,r) = (10.^(xi(D(cond).std)));
                end
            end
            if options.llh.original
                if cond == 11 || cond == 12
                    sigma2(1,6) = sigma2(1,6) + D(cond).u(3)*10.^(2*xi(113));
                end
            end
            if nargout > 1
                dsigma2 = zeros(1,20,numel(xi),size(D(cond).my,3));
                for iobs = 1:20
                    dsigma2(1,iobs,D(cond).std(iobs),:) = dsigma2(1,iobs,D(cond).std(iobs),:) + ...
                        10.^(xi(D(cond).std(iobs)))*log(10);
                end
            end
            switch options.llh.distribution
                case 'log-laplace' %Note: factor -log(D(cond).my) neglected
                    error('to do')
                    %                     % check gradient!
                    %                     y_ch = bsxfun(@minus,log(D(cond).my),log(sol(cond).y));
                    %                     nlogL = nlogL + sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),log(2*sigma2))+...
                    %                         bsxfun(@rdivide,abs(y_ch),sigma2),1),3),2);
                    %                     if nargout > 1
                    %                         dnlogL = dnlogL + permute(sum(sum(nansum(...
                    %                             bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,abs(y_ch),sigma2)),...
                    %                             sigma2),[1,2,4,3]),dsigma2)-...
                    %                             bsxfun(@times,permute(bsxfun(@rdivide,sign(y_ch),sigma2),[1,2,4,3]),...
                    %                             bsxfun(@rdivide,sol(cond).sy,sol(cond).y))...
                    %                             ,1),4),2),[3,2,1]);
                    %                     end
                case 'log10-normal'  %Note: factor -log(D(cond).my) neglected
%                     if cond == 2
%                         my = D(cond).my;
%                         my(:,[1:10,12:end],:) = log10(my(:,[1:10,12:end],:));
%                         y = sol(cond).y;
%                         y(:,[1:10,12:end]) = log10(y(:,[1:10,12:end]));
%                         y_ch = bsxfun(@minus,my,y);
%                     else
                        y_ch = bsxfun(@minus,log10(D(cond).my),log10(sol(cond).y));
                   % end
                    nlogL = nlogL + 0.5*(sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),...
                        log(2*pi*sigma2))+...
                        bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2),1),3),2));
                    if nargout > 1
                        dnlogL = dnlogL + 0.5*(permute(sum(sum(nansum(...
                            bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                            bsxfun(@power,y_ch,2),sigma2)),sigma2),[1,2,4,3]),dsigma2) -...
                            bsxfun(@times,1/log(10)*permute(2*bsxfun(@rdivide,y_ch,sigma2),[1,2,4,3]),...
                            bsxfun(@rdivide,sol(cond).sy,sol(cond).y))...
                            ,1),4),2),[3,2,1]));
                    end
                case 'log-normal'  %Note: factor -log(D(cond).my) neglected
                    y_ch = bsxfun(@minus,log(D(cond).my),log(sol(cond).y));
                    if options.llh.lsqnonlin
                        for ir = 1:size(D(cond).my,3)
                            temp = bsxfun(@rdivide,y_ch(:,:,ir),sqrt(sigma2(:,:,ir)));
                            res = [res;temp(:)];
                            temp_res_err = sqrt(bsxfun(@times,~isnan(D(cond).my(:,:,ir)),log(2*pi*sigma2(:,:,ir))+50));
                            res_err = temp_res_err(:);
                            ind = find(res_err>0);
                            res = [res;res_err(ind)];
                            if nargout > 1
                                stemp = bsxfun(@minus,-bsxfun(@times,temp./sqrt(sigma2(:,:,ir)),...
                                    bsxfun(@rdivide,dsigma2(:,:,:,ir),2*sqrt(sigma2(:,:,ir)))),...
                                    bsxfun(@times,1./(bsxfun(@times,sqrt(sigma2(:,:,ir)),sol(cond).y)),sol(cond).sy));
                                stemp = reshape(stemp,numel(temp),numel(xi));
                                sres = [sres; stemp];
                                
                                temp_sres_err = bsxfun(@times,bsxfun(@rdivide,1./temp_res_err,2*sigma2(:,:,ir)),...
                                    dsigma2(:,:,:,ir));
                                temp_sres_err = reshape( temp_sres_err,numel(temp_res_err),numel(xi));
                                sres = [sres; temp_sres_err(ind,:)];
                            end
                        end
                    else
                        nlogL = nlogL + 0.5*(sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),log(2*pi*sigma2))+...
                            bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2),1),3),2));
                        %  nlogL = nlogL + 0.5*(sum(sum(nansum(...
                        %    bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2),1),3),2));
                        if nargout > 1
                            dnlogL = dnlogL + 0.5*(permute(sum(sum(nansum(...
                                bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                                bsxfun(@power,y_ch,2),sigma2)),sigma2),[1,2,4,3]),dsigma2) -...
                                bsxfun(@times,permute(2*bsxfun(@rdivide,y_ch,sigma2),[1,2,4,3]),...
                                bsxfun(@rdivide,sol(cond).sy,sol(cond).y))...
                                ,1),4),2),[3,2,1]));
                        end
                    end
                case 'laplace'
                    y_ch = bsxfun(@minus,D(cond).my,sol(cond).y);
                    nlogL = nlogL + sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),log(2*sigma2))+...
                        bsxfun(@rdivide,abs(y_ch),sigma2),1),3),2);
                    if nargout > 1
                        dnlogL = dnlogL + permute(sum(sum(nansum(...
                            bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,abs(y_ch),sigma2)),...
                            sigma2),[1,2,4,3]),dsigma2)-...
                            bsxfun(@times,permute(bsxfun(@rdivide,sign(y_ch),sigma2),[1,2,4,3]),sol(cond).sy)...
                            ,1),4),2),[3,2,1]);
                    end
                case 'normal'
                    y_ch = bsxfun(@minus,D(cond).my,sol(cond).y);
                    nlogL = nlogL + 0.5*(sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),log(2*pi*sigma2))+...
                        bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2),1),3),2));
                    if nargout > 1
                        dnlogL = dnlogL + 0.5*(permute(sum(sum(nansum(...
                            bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                            bsxfun(@power,y_ch,2),sigma2)),sigma2),[1,2,4,3]),dsigma2) -...
                            bsxfun(@times,permute(2*bsxfun(@rdivide,y_ch,sigma2),[1,2,4,3]),sol(cond).sy)...
                            ,1),4),2),[3,2,1]));
                    end
            end
        end
end

if isfield(options.llh,'parameter_prior')
    if options.llh.lsqnonlin
        error('todo: include res for prior')
    else
        nlogL = nlogL + 0.5*nansum((xi-options.llh.parameter_prior.mean).^2./ options.llh.parameter_prior.sigma2);
        if nargout > 1
            dnlogL = nansum([dnlogL,(xi-options.llh.parameter_prior.mean)./ options.llh.parameter_prior.sigma2],2);
            
        end
    end
end

if options.llh.lsqnonlin
    ind = find(~isnan(res));
    varargout{1} = res(ind);
    if nderiv>=1
        if options.llh.reduced_woinit
            varargout{2} = sres(ind,[1:24,27:end]);
        else
            varargout{2} = sres(ind,:);
        end
    end
else
    varargout{1} = nlogL;
    if nderiv>=1
        if options.llh.reduced_woinit
            varargout{2} =dnlogL([1:24,27:end]);% squeeze(dsigma2(1,1,[1:24,27:end]));%
        else
            varargout{2} = dnlogL;
        end
    end
end