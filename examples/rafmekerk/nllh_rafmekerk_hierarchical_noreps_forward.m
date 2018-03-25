function [varargout] = nllh_rafmekerk_hierarchical_noreps_forward(theta,D,scOptions)

amiOptions = amioption();
simfun = @simulate_rafmekerk_hierarchical_noreps_forward;

switch nargout
    case 1
        [varargout{1}] = hieropt_nllh_forward(simfun,theta,D,amiOptions,scOptions);
    case 2
        [varargout{1},varargout{2}] = hieropt_nllh_forward(simfun,theta,D,amiOptions,scOptions);
    case 3
        [varargout{1},varargout{2},varargout{3}] = hieropt_nllh_forward(simfun,theta,D,amiOptions,scOptions);
end

end