function [ varargout ] = hieropt_nllh_forward( varargin )
% hieropt_nllh_foroward uses the hierarchical approach to compute the nllh.
% The first derivative is computed using the forward approach, additionally
% a FIM approximation to the second derivative can be obtained as third
% return value.
% The function can be run in two modes: Either a simulation function is
% passed (b_sim = true), or only the simulated values and the optimal
% scalings (obtained from hieropt_scalings) are passed.
% For an example usage see
% examples/HierOpt_Examples/jakstat_small/nllh_jakstat_hierarchical
% (_offsets).
%
% Parameters:
%   b_sim
%   if b_sim == true
%     simfun
%     theta
%     D
%     amiOptions
%     scOptions
%   else
%     sim
%     D
%     b
%     c
%     sigma2
%
% Return Values:
%  [nllh,snllh,s2nllh]
%
% History:
%   2018/01/12: Yannik Schaelte

% simulate first?
b_sim = varargin{1};

if b_sim
    simfun = varargin{2};
    theta = varargin{3};
    D = varargin{4};
    amiOptions = varargin{5};
    scOptions = varargin{6};

    switch nargout
        case 1
            [varargout{1}] = hieropt_nllh_forward_withsim(simfun,theta,D,amiOptions,scOptions);
        case 2
            [varargout{1},varargout{2}] = hieropt_nllh_forward_withsim(simfun,theta,D,amiOptions,scOptions);
        case 3
            [varargout{1},varargout{2},varargout{3}] = hieropt_nllh_forward_withsim(simfun,theta,D,amiOptions,scOptions);
    end
else
    sim = varargin{2};
    D = varargin{3};
    b = varargin{4};
    c = varargin{5};
    sigma2 = varargin{6};
    
    switch nargout
        case 1
            [varargout{1}] = hieropt_nllh_forward_withoutsim(sim,D,b,c,sigma2);
        case 2
            [varargout{1},varargout{2}] = hieropt_nllh_forward_withoutsim(sim,D,b,c,sigma2);
        case 3
            [varargout{1},varargout{2},varargout{3}] = hieropt_nllh_forward_withoutsim(sim,D,b,c,sigma2);
    end
end

end % function

function [ varargout ] = hieropt_nllh_forward_withsim(simfun,theta,D,amiOptions,scOptions)

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
    sol = simfun(D(ie).t,theta,D(ie).k,[],amiOptions);
    
    if (sol.status ~= 0)
        error('hieropt: could not integrate ODE.');
    end
    
    sim(ie).y = sol.y;
    if nargout > 1
        sim(ie).sy = sol.sy;
    end
end

% optimal scalings
[b,c,sigma2] = hieropt_scalings(sim,D,scOptions);

% nllh, grad, fim computed from given data
switch nargout
    case 1
        [varargout{1}] = hieropt_nllh_forward_withoutsim(sim,D,b,c,sigma2);
    case 2
        [varargout{1},varargout{2}] = hieropt_nllh_forward_withoutsim(sim,D,b,c,sigma2);
    case 3
        [varargout{1},varargout{2},varargout{3}] = hieropt_nllh_forward_withoutsim(sim,D,b,c,sigma2);
end

end % function

function [ varargout ] = hieropt_nllh_forward_withoutsim(sim,D,b,c,sigma2)

ne = size(D,2);
if nargout > 1
    n_theta = size(sim(1).sy,3);
end

% initialize output
nllh = 0;
if nargout > 1
    grad = zeros(n_theta,1);
    if nargout > 2
        fim = zeros(n_theta);
    end
end

for ie = 1:ne
    nr = size(D(ie).Y,3);
    
    b_e = b{ie};
    c_e = c{ie};
    sigma2_e = sigma2{ie};
    
    y_ch = bsxfun(@minus,D(ie).Y,bsxfun(@times,c_e,sim(ie).y)+b_e);
    
    nllh = nllh + 0.5*sum(sum(nansum(bsxfun(@times,~isnan(D(ie).Y),log(2*pi*sigma2_e))+...
        bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2_e),1),3),2);
    
    if nargout > 1
        dy_ch = - bsxfun(@times,permute(c_e,[1,2,4,3]),repmat(sim(ie).sy,[1,1,1,nr]));
        
        grad = grad + permute(sum(sum(nansum(...
            bsxfun(@times,permute(...
            bsxfun(@rdivide,y_ch,sigma2_e),[1,2,4,3]),...
            dy_ch),1),4),2),[3,2,1]);
        
        if nargout > 2
            for j=1:n_theta
                for k=1:n_theta
                    fim(j,k) = fim(j,k) + nansum(nansum(nansum(...
                        bsxfun(@rdivide,...
                        bsxfun(@times,dy_ch(:,:,:,j),dy_ch(:,:,:,k)),...
                        sigma2_e))));
                end
            end
        end
    end
end

varargout{1} = nllh;
if nargout > 1
    varargout{2} = grad;
    if nargout > 2
        varargout{3} = fim;
    end
end

end % function
