function [varargout] = nllh_rafmekerk_standard(theta, D, scOptions)

switch scOptions.distribution
    case 'normal'
        simfun = @simulate_rafmekerk_standard;
    case 'laplace'
        simfun = @simulate_rafmekerk_laplace_standard;
end
    
    %% Objective Function Evaluation
    
    % Initialization
    llh = 0;
    sllh = zeros(28, 1);
    
    % Integration
    if (nargout == 1)
        amiOptions.sensi = 0;
        for j = 1 : 3
            sol = simfun(D(j).t, theta, D(j).condition, D(j), amiOptions);
            if sol.status < 0
                llh = -inf;
            end
            llh = llh + sol.llh;
        end
        
    elseif (nargout == 2)
        amiOptions.sensi = 1;
        for j = 1 : 3
            sol = simfun(D(j).t, theta, D(j).condition, D(j), amiOptions);
            if sol.status < 0
                llh = -inf;
            end
            llh = llh + sol.llh;
            sllh = sllh + sol.sllh;
        end
        
    end
    
    % Assignment
    switch (nargout)
        case{0,1}
            varargout{1} = -llh;
            
        case 2
            varargout{1} = -llh;
            varargout{2} = -sllh;

    end

end