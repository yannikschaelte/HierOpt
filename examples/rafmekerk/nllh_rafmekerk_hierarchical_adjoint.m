function [varargout] = nllh_rafmekerk_hierarchical_adjoint(theta,D,scOptions)

amiOptions = amioption();
switch scOptions.distribution
    case 'normal'
        simfun = @simulate_rafmekerk_hierarchical_adjoint;
    case 'laplace'
        simfun = @simulate_rafmekerk_laplace_hierarchical_adjoint;
end

switch nargout
    case 1
        [varargout{1}] = hieropt_nllh_adjoint(simfun,theta,D,amiOptions,scOptions);
    case 2
        [varargout{1},varargout{2}] = hieropt_nllh_adjoint(simfun,theta,D,amiOptions,scOptions);
    case 3
        [varargout{1},varargout{2},varargout{3}] = hieropt_nllh_adjoint(simfun,theta,D,amiOptions,scOptions);
end

end
