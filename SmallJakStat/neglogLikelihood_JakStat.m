function [varargout] = neglogLikelihood_JakStat(varargin)
% compute the loglikelihood, gradient and Hessian for JakStat data

%INPUT:
%xi: parameter vector
%kappa: konstants vector
%D: data struct of the form D.t D.my
%distr: indicates distribution of measurement noise,
%           * normal
%           * laplace
%save_pv: indicates weather the proportionality factors and variances are
%         saved in a .mat file (1) or not (0)
%initpar: indicates if init_STAT is a parameter or fixed
%           * 0: parameter
%           * 1: fixed


%OUTPUT
%negative-log-Likelihood
%Gradient of the negative-log-Likelihood
%Approximation of the Hessian of the neg-log-Likelihood (FIM)


% CHECK/ASSIGN INPUTS:
if nargin >= 3
    xi = varargin{1};
    kappa = varargin{2};
    D = varargin{3};
    distr = varargin{4};
    options = varargin{5};
else
    error('Not enough input arguments!');
end

try
    save_pv = options.llh.save_analytical;
    options.ami.atol = 1e-12;
    options.ami.rtol = 1e-12;
    if nargout>1
        options.ami.sensi = 1;
    else
        options.ami.sensi = 0;
    end
    
    %% SIMULATION
    switch options.llh.approach
        case 'hierarchical'
            sol = simulate_JakStat_hierarchical(D.t,xi(1:11),kappa,[],options.ami);
            simulation(1).y = sol.y;
            if nargout > 1
                simulation(1).sy = sol.sy;
            end
            %% NEGATIVE-LOG-LIKELIHOOD, GRADIENT, FIM
            if nargout > 2
                [nlLH, gradnlLH, FIMnlLH] = nlLH_fgh(simulation,D,distr,options.llh,save_pv);
            elseif nargout > 1
                [nlLH, gradnlLH] = nlLH_fgh(simulation,D,distr,options.llh,save_pv);
            else
                nlLH = nlLH_fgh(simulation,D,distr,options.llh,save_pv);
            end
        case 'standard'
            sol = simulate_JakStat(D.t,xi,kappa,[],options.ami);
            my = sol.y;
            if nargout>1
                dmydxi = sol.sy;
            end
            %% LOGLIKELIHOOD, GRAD
            % initialize logL, grad,fish
            logL = 0;
            if nargout>1
                grad = zeros(length(xi),1);
                if nargout>2
                    fish = zeros(length(xi),length(xi));
                end
            end
            nt = size(my,1); %number of timepoins
            no = size(my,2); %number of outputs
            np = length(xi)-no; % number of dynamic parameters
            
            sigma = zeros(no,1);
            for i = 1:no
                sigma(i) = 10^xi(np+i);
            end
            
            D.sigma_my = repmat(sigma',[nt,1]);
            resmy = reshape((my-D.my)./(D.sigma_my),nt*no,1);
            
            switch distr
                case 'normal'
                    if nargout>1
                        sresmy = reshape(bsxfun(@times,dmydxi,1./(D.sigma_my)),nt*no,np);
                        resmyc = reshape((my-D.my),nt*no,1);
                    end
                    for i = 1:no
                        logL = logL - 0.5*nansum(log(2*pi*(sigma(i)).^2) + resmy((i-1)*nt+1:i*nt).^2);
                    end
                    if nargout>1
                        % Compute Gradient with sensitivities
                        dsigmadxi = zeros(length(xi),no);
                        
                        for i = 1:no
                            dsigmadxi(np+i,i) = sigma(i)*log(10);
                            grad(1:np) = grad(1:np) - nansum(bsxfun(@times,dmydxi(:,i+(0:no:np*no-i)),((resmyc((i-1)*nt+1:i*nt)))./(sigma(i).^2)))';
                            grad = grad - nansum((1-resmy((i-1)*nt+1:i*nt).^2))*1/(sigma(i))*dsigmadxi(:,i);
                            
                        end
                        if nargout>2
                            for i=1:no
                                dsigma2dxi = 2*sigma(i)*dsigmadxi(:,i);
                                d2sigma2dxi2 = 2*(dsigmadxi(:,i)*dsigmadxi(:,i)');
                                d2sigma2dxi2(np+i,np+i) = d2sigma2dxi2(np+i,np+i) + 2*sigma(i)^2*log(10)^2;
                                G = resmy((i-1)*nt+1:i*nt).^2;
                                Happ9_1 = nansum((1/(2*sigma(i)^4))*(1-2*G))*(dsigma2dxi*dsigma2dxi');
                                Happ9_2 = -1*nansum((1/(2*sigma(i)^2))*(1-G))*d2sigma2dxi2;
                                Happ78 = -1*dsigma2dxi*nansum(bsxfun(@times,[-1*dmydxi(:,i+(0:no:np*no-i)),zeros(nt,no)],((resmyc((i-1)*nt+1:i*nt)))))./(sigma(i).^4);
                                Happ36 = Happ78';
                                Happ1245 = -1*[-1*dmydxi(:,i+(0:no:np*no-i)),zeros(nt,no)]'*[-1*dmydxi(:,i+(0:no:np*no-i)),zeros(nt,no)]*1/(sigma(i)^2);
                                fish = fish + Happ9_1+Happ9_2+Happ78+Happ36+Happ1245;
                            end
                        end
                        
                    end
                case 'laplace'
                    resmyc = reshape((my-D.my),nt*no,1);
                    for i = 1:no
                        logL = logL - nansum(log(2*sigma(i)) + abs(resmyc((i-1)*nt+1:i*nt))./sigma(i));
                    end
                    if nargout>1
                        dbdxi = zeros(length(xi),no);
                        
                        for i = 1:no
                            dbdxi(np+i,i) = sigma(i)*log(10);
                            grad(1:np) = grad(1:np) - nansum(bsxfun(@times,(1/sigma(i))*dmydxi(:,i+(0:no:np*no-i)),sign(resmyc((i-1)*nt+1:i*nt))))';
                            grad = grad + nansum(-(1/(sigma(i)))+ abs(resmyc((i-1)*nt+1:i*nt))./(sigma(i)^2))*dbdxi(:,i);
                        end
                        if nargout>2
                            for i = 1:no
                                d2bdxi2 = zeros(length(xi),length(xi));
                                d2bdxi2(np+i,np+i) = sigma(i)*(log(10)^2);
                                Happ91 = nansum(-(1/sigma(i))+ abs(resmyc((i-1)*nt+1:i*nt))./(sigma(i)^2))*d2bdxi2;
                                Happ92 = nansum((1/sigma(i)^2)-2*abs(resmyc((i-1)*nt+1:i*nt))./(sigma(i)^3))*dbdxi(:,i)*dbdxi(:,i)';
                                Happ78 = zeros(length(xi),length(xi));
                                Happ78(1:length(xi),1:np) = dbdxi(:,i)*nansum(-1*bsxfun(@times,dmydxi(:,i+(0:no:np*no-i)),sign(resmyc((i-1)*nt+1:i*nt))./(sigma(i)^2)));
                                Happ36 = Happ78';
                                fish = fish +Happ91+Happ92-Happ78-Happ36;
                            end
                            
                        end
                    end
            end
    end
catch error_thrown
    warning(['Evaluation of likelihood failed. ',error_thrown.message]);
    nlLH = nan;
    gradnlLH = nan(length(xi),1);
    FIMnlLH = nan(length(xi),length(xi));
    logL = nan;
    grad = nan(length(xi),1);
    fish = nan(length(xi),length(xi));
end
if sol.status < 0
    warning('Failed to integrate ODE.')
    nlLH = nan;
    gradnlLH = nan(length(xi),1);
    FIMnlLH = nan(length(xi),length(xi));
    logL = nan;
    grad = nan(length(xi),1);
    fish = nan(length(xi),length(xi));
end
switch options.llh.approach
    case 'hierarchical'
        switch nargout
            case{0,1}
                varargout{1} = nlLH;
            case 2
                varargout{1} = nlLH;
                varargout{2} = gradnlLH;
            case 3
                varargout{1} =  nlLH;
                varargout{2} =  gradnlLH;
                varargout{3} =  FIMnlLH;
        end
    case 'standard'
        switch nargout
            case{0,1}
                varargout{1} = -logL;
            case 2
                varargout{1} = -logL;
                varargout{2} = -grad;
            case 3
                varargout{1} =  -logL;
                varargout{2} =  -grad;
                varargout{3} =  -fish;
        end
end
end

