function [ varargout ]  = nllh_jakstat_hierarchical(theta,D,scOptions)

amiOptions.rtol = 1e-10;
amiOptions.atol = 1e-10;
simfun = @simulate_jakstat_hierarchical;

switch nargout
    case 1
        [varargout{1}] = hieropt_nllh_forward(true,simfun,...
            theta,D,amiOptions,scOptions);
    case 2
        [varargout{1},varargout{2}] = hieropt_nllh_forward(true,simfun,...
            theta,D,amiOptions,scOptions);
    case 3
        [varargout{1},varargout{2},varargout{3}] = hieropt_nllh_forward(true,simfun,...
            @simulate_jakstat_hierarchical,...
            theta,D,amiOptions,scOptions);
end

end

% amiOptions.sensi_meth = 'forward';
% if nargout == 1
%     amiOptions.sensi = 0;
% else
%     amiOptions.sensi = 1;
% end
%
% % for every experiment and replicate, do the optimization, then aggregate
% n_e = size(D,2);
%
% % forward simulation
%
% sim = struct([]);
% for ie = 1:n_e
%     sol = simulate_jakstat_hierarchical(D(ie).t,theta,kappa(:,ie),[],amiOptions);
%
%     if (sol.status ~= 0)
%         error('Could not integrate ODE.');
%     end
%
%     sim(ie).y = sol.y;
%     if nargout > 1
%         sim(ie).sy = sol.sy;
%     end
% end
%
% [b,c,sigma2] = hieropt_scalings(sim,D,scOptions);
%
% switch nargout
%     case 1
% %         nllh = nlLH_fgh(sim,D,'normal',scOptions,false);
%         nllh = opt_nllh_forward(false,sim,D,b,c,sigma2);
%         varargout{1} = nllh;
%     case 2
% %         [nllh,grad] = nlLH_fgh(sim,D,'normal',scOptions,false);
%         [nllh,grad] = opt_nllh_forward(false,sim,D,b,c,sigma2);
%         varargout{1} = nllh;
%         varargout{2} = grad;
%     case 3
% %         [nllh,grad,fim] = nlLH_fgh(sim,D,'normal',scOptions,false);
%         [nllh,grad,fim] = opt_nllh_forward(false,sim,D,b,c,sigma2);
%         varargout{1} = nllh;
%         varargout{2} = grad;
%         varargout{3} = fim;
% end
%
% end
