 function [varargout] = neglogLikelihood_RafMekErk_standard(varargin)
% compute the negative loglikelihood and gradient for RafMekErk data

%INPUT:
%xi: parameter vector
%D: data struct of the form D.t D.my
%distr: indicates distribution of measurement noise,
%           * normal
%           * laplace

%OUTPUT
%loglikelihood
%gradient

% CHECK/ASSIGN INPUTS:
if nargin >= 3
    xi = varargin{1};
    D = varargin{2};
    distr = varargin{3};
else
    error('Not enough input arguments!');
end
n_theta = 20; %number of dynamical parameters

if nargin == 4
    prior2 = varargin{4};
end

try
    %global errorCount
    u = D.conditions;
    n_u = size(u,1);
    n_sigma = length(xi)-n_theta;
    n_r = size(D.measurement{1},2);    
    theta = xi(1:n_theta);% (log10) parameters
    noise2 = bsxfun(@times,ones(length(D.t{1}),n_r),xi(n_theta+1:n_theta+n_sigma)');    
    J = 0;
    if nargout>1
        dJdtheta = zeros(length(xi),1);
        options_simu.sensi = 1;
    else
        options_simu.sensi = 0;
    end

    for i = 1:n_u
        t_m = D.t{i};
        n_t = length(t_m);
        sol = simulate_RafMekErk(t_m,theta,u(i,:),[],options_simu);
        %[status,~,~,y_i,~,sy_i] = simulate_RafMekErklin_ana(t_m,theta,u(i,:));

        if sol.status < 0
            error(['failed to integrate ODE for experiment ' num2str(i)])
        end

        y = sol.y(:);
        S_theta = reshape(sol.sy,[],n_theta);
        sigma2 = bsxfun(@power,10,noise2);
        if nargout > 1
            dsigma2dtheta = zeros(length(y),n_sigma);           
            dsigma2dtheta(1:n_t,1) = log(10)*ones(n_t,1);
            dsigma2dtheta(n_t+1:2*n_t,2) = log(10)*ones(n_t,1);
            dsigma2dtheta(2*n_t+1:3*n_t,3) = log(10)*ones(n_t,1);
            dsigma2dtheta(3*n_t+1:4*n_t,4) = log(10)*ones(n_t,1);
            dsigma2dtheta(4*n_t+1:5*n_t,5) = log(10)*ones(n_t,1);
            dsigma2dtheta(5*n_t+1:6*n_t,6) = log(10)*ones(n_t,1);
            dsigma2dtheta(6*n_t+1:7*n_t,7) = log(10)*ones(n_t,1);
            dsigma2dtheta(7*n_t+1:8*n_t,8) = log(10)*ones(n_t,1);

            dsigmadtheta = [zeros(length(y),n_theta-n_r),...
                zeros(length(y),n_r),bsxfun(@times,dsigma2dtheta,sigma2(:))];

            % Parameter derivative
            dyxi_ik = S_theta;
            dyxi_ik_t = [dyxi_ik, zeros(length(y),n_sigma)];
        end

        %% Objective function (J)
        if nargout>1
            switch distr
                case 'normal'
                    [J_i,grad_J] = J_normal(D.measurement{i}(:),y,sigma2(:),dyxi_ik_t,dsigmadtheta);
                case 'laplace'
                    [J_i,grad_J] = J_laplace(D.measurement{i}(:),y,sigma2(:),dyxi_ik_t,dsigmadtheta);
            end
            J = J + J_i;
            dJdtheta = dJdtheta + grad_J;
        else
            switch distr
                case 'normal'
                    J_i = J_normal(D.measurement{i}(:),y,sigma2(:),[],[]);
                case 'laplace'
                    J_i = J_laplace(D.measurement{i}(:),y,sigma2(:),[],[]);
            end
            J = J + J_i;
        end
    end
    
    if(nargin == 4)
        J = J + 0.5*sum(xi(1:12).^2)/prior2;
        if nargout > 1
            dJdtheta = dJdtheta + [xi(1:12)/prior2;zeros(16,1)];
        end
    end

catch error_thrown
    warning(['Evaluation of likelihood failed. ',error_thrown.message]);
    J = inf;
    dJdtheta = zeros(length(xi),1);
end
switch nargout
    case{0,1}
        varargout{1} = J;
    case 2
        varargout{1} = J;
        varargout{2} = dJdtheta;
end

 end

 
 function [varargout] = J_laplace(y_m,y,sigma2,dydtheta,dsigma2dtheta)
%OUTPUT
%loglikelihood: J
%gradient: dJdtheta

% calculation of loglikelihood value for non NaNs 
ind = find(~isnan(y_m));
J=sum((log(2*sigma2(ind))+abs(y_m(ind)-y(ind))./sigma2(ind)));

if nargout > 1
    %calculation of gradient with respect to theta
    grad_J = zeros(size(y_m));
    grad_J(ind) = -sign(y_m(ind)-y(ind))./sigma2(ind);

    % calculation of gradient with respect to sigma^2
    grad_J_sig = zeros(size(y_m));
    grad_J_sig(ind) = 1./sigma2(ind).*(1-abs(y_m(ind)-y(ind))./sigma2(ind));

    % grad_J = grad_J*dydtheta + grad_J*dydsigma2*dsigma2dtheta
    dJdtheta = (grad_J' * dydtheta + grad_J_sig' * dsigma2dtheta)';
end

switch nargout
    case{0,1}
        varargout{1} = J;
    case 2
        varargout{1} = J;
        varargout{2} = dJdtheta;
end

 end

 function [varargout] = J_normal(y_m,y,sigma2,dydtheta,dsigma2dtheta)
%OUTPUT
%loglikelihood: J
%gradient: dJdtheta

% calculation of loglikelihood value for non NaNs 
ind = find(~isnan(y_m));
J=sum(0.5*(log(2*pi*sigma2(ind))+(y_m(ind)-y(ind)).^2./sigma2(ind)));

if nargout > 1
    %calculation of gradient with respect to theta
    grad_J = zeros(size(y_m));
    grad_J(ind) = -(y_m(ind)-y(ind))./sigma2(ind);

    % calculation of gradient with respect to sigma^2
    grad_J_sig = zeros(size(y_m));
    grad_J_sig(ind) = 0.5*1./sigma2(ind).*(1-(y_m(ind)-y(ind)).^2./sigma2(ind));

    % grad_J = grad_J*dydtheta + grad_J*dydsigma2*dsigma2dtheta
    dJdtheta = (grad_J' * dydtheta + grad_J_sig' * dsigma2dtheta)';
end

switch nargout
    case{0,1}
        varargout{1} = J;
    case 2
        varargout{1} = J;
        varargout{2} = dJdtheta;
end

end