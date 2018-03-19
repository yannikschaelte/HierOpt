function [ varargout ] = nlllh_erbb_standard(theta, D)

amiOptions.rtol = 1e-8;
amiOptions.atol = 1e-16;
amiOptions.maxsteps = 2e5;
amiOptions.interpType = 2;
amiOptions.sensi_meth = 'adjoint';

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
        amiOptions.sensi = 2;
    otherwise
        error('Only supports up to 3 outputs.');
end

for ie = 1:n_e
    n_r = size(D(ie).Y,3);
    for ir = 1:n_r
        % set data for simulation
        amiData.t = D(ie).t;
        amiData.Y = D(ie).Y(:,:,ir);
        amiData.Sigma_Y = D(ie).Sigma_Y(:,:,ir);
        amiData.condition = D(ie).condition(:);
        amiData = amidata(amiData);
        
        sol = simulate_erbb_standard([], theta, [], amiData, amiOptions);
        
        if sol.status ~= 0
            error('Could not integrate ODE.');
        end
        
        varargout{1} = varargout{1} - sol.llh;
        if nargout > 1
            varargout{2} = varargout{2} - sol.sllh;
            if nargout > 2
                varargout{3} = varargout{3} - sol.s2llh;
            end
        end
    end 
end

end

