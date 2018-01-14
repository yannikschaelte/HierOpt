function [varargout] = nllh_rafmekerk_hierarchical(theta,D,scOptions)

amiOptions = amioption();
simfun = @simulate_rafmekerk_hierarchical;

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