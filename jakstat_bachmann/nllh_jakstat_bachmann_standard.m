function [ varargout ] = nllh_jakstat_bachmann_standard(theta, D)

amiOptions.sensi_meth = 'forward';

switch nargout
    case 1
        nllh = 0;
        amiOptions.sensi = 0;
    case 2
        nllh = 0;
        snllh = 0;
        amiOptions.sensi = 1;
    case 3
        nllh = 0;
        snllh = 0;
        s2nllh = 0;
        amiOptions.sensi = 2;
    otherwise
        error('Only supports up to 3 outputs.');
end

% for every experiment and replicate, do the optimization
n_e = size(D,2);

for ie = 1:n_e
    n_r = size(D(ie).Y,3);
    for ir = 1:n_r
        clear('amiData');
        amiData.t = D(ie).t;
        amiData.Y = D(ie).Y(:,:,ir);
        amiData.condition = D(ie).k(:);
        amiData = amidata(amiData);
        
        sol = simulate_jakstat_bachmann([], theta, [], amiData, amiOptions);
    end
end

end

