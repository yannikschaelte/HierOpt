function [ varargout ]  = nllh_jakstat_hierarchical(theta,kappa,D,scOptions)

amiOptions.rtol = 1e-10;
amiOptions.atol = 1e-10;
amiOptions.sensi_meth = 'forward';

% for every experiment and replicate, do the optimization
n_e = size(D,2);

% prepare output and sensi
switch nargout
    case 1
        varargout{1} = 0;
        amiOptions.sensi = 0;
    case 2
        varargout{1} = 0;
        varargout{2} = 0;
        amiOptions.sensi = 1;
    case 3
        varargout{1} = 0;
        varargout{2} = 0;
        varargout{3} = 0;
        amiOptions.sensi = 1;
    otherwise
        error('Only supports up to 3 outputs.');
end

% forward simulation

sim = struct([]);
for ie = 1:n_e
    
    sol = simulate_jakstat_hierarchical(D(ie).t,theta,kappa(:,ie),[],amiOptions);
    
    if (sol.status ~= 0)
        error('Could not integrate ODE.');
    end

    sim(ie).y = sol.y;
    if nargout > 1
        sim(ie).sy = sol.sy;
    end
end

if nargout == 1
    nllh = nlLH_fgh(sim,D,'normal',scOptions,false);
    varargout{1} = nllh;
else
    [nllh,snllh] = nlLH_fgh(sim,D,'normal',scOptions,false);
    varargout{1} = nllh;
    varargout{2} = snllh;
end

end
