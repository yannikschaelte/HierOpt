function [ varargout ] = nllh_erbb_hierarchical_forward(theta, D)

amiOptions.rtol = 1e-8;
amiOptions.atol = 1e-16;
amiOptions.maxsteps = 2e5;
amiOptions.interpType = 2;
simfun = @simulate_erbb_hierarchical_forward;

switch nargout
    case 1
        [varargout{1}] = hieropt_nllh_forward(simfun,theta,D,amiOptions,scOptions);
    case 2
        [varargout{1},varargout{2}] = hieropt_nllh_forward(simfun,theta,D,amiOptions,scOptions);
    case 3
        [varargout{1},varargout{2},varargout{3}] = hieropt_nllh_forward(simfun,theta,D,amiOptions,scOptions);
end

end