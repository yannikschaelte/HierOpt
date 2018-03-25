function [ varargout ] = hieropt_nllh_forward( simfun, theta, D, amiOptions, scOptions )
% Uses the hierarchical approach to compute the nllh.
%
% The first derivative is computed using the forward approach, additionally
% a FIM approximation to the second derivative can be obtained as third
% return value.
% The function can be run in two modes: Either a simulation function is
% passed (b_sim = true), or only the simulated values and the optimal
% scalings (obtained from hieropt_scalings) are passed.
%
% The objective function interface must be that of the DE integration 
% toolbox AMICI.
%
% Input:
%   simfun
%   theta
%   D
%   amiOptions
%   scOptions
%
% Output:
%  [nllh,snllh,s2nllh]
%
% History:
%   2018/01/12: Yannik Schaelte

% simulate

% set correct amiOptions
amiOptions.sensi_meth = 'forward';
if nargout == 1
    amiOptions.sensi = 0;
else
    amiOptions.sensi = 1;
end

ne = size(D,2);

% forward simulation
sim = struct([]);
for ie = 1:ne
    sol = simfun(D(ie).t,theta,D(ie).condition,[],amiOptions);
    
    if (sol.status ~= 0)
        error('hieropt:nllh', "Could not integrate ODE.");
    end
    
    sim(ie).y = sol.y;
    if nargout > 1
        sim(ie).sy = sol.sy;
    end
end

switch nargout
    case 1
        [varargout{1}] = hieropt_nllh_nosim(sim,D,scOptions);
    case 2
        [varargout{1},varargout{2}] = hieropt_nllh_nosim(sim,D,scOptions);
    case 3
        [varargout{1},varargout{2},varargout{3}] = hieropt_nllh_nosim(sim,D,scOptions);
end

end % function