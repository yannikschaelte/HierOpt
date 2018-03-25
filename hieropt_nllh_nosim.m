function [ varargout ] = hieropt_nllh_nosim( sim, D, varargin )
% Input:
%   sim
%   D
%   <either>
%     scOptions
%   <or>
%     b
%     c
%     noise
%     distribution
%
% Output:

if nargin == 3
    scOptions = varargin{1};
    [b,c,noise] = hieropt_scalings(sim,D,scOptions);
    distribution = scOptions.distribution;
else
    b = varargin{1};
    c = varargin{2};
    noise = varargin{3};
    distribution = varargin{4};
end

switch distribution
    case 'normal'
        switch nargout
            case 1
                [varargout{1}] = nllh_normal(sim,D,b,c,noise);
            case 2
                [varargout{1},varargout{2}] = nllh_normal(sim,D,b,c,noise);
            case 3
                [varargout{1},varargout{2},varargout{3}] = nllh_normal(sim,D,b,c,noise);
        end
    case 'laplace'
        switch nargout
            case 1
                [varargout{1}] = nllh_laplace(sim,D,b,c,noise);
            case 2
                [varargout{1}] = nllh_laplace(sim,D,b,c,noise);
            case 3
                error('hieropt:nllh_nosim',"FIM not supported for Laplace errors");
        end
end

end


function [ varargout ] = nllh_normal( sim, D, b, c, noise )

n_e = size(D,2);
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

% compute output
for ie = 1:n_e
    n_r = size(D(ie).Y,3);
    b_e = b{ie};
    c_e = c{ie};
    noise_e = noise{ie};
    
    y_e = D(ie).Y;
    h_e = repmat(sim(ie).y, [1,1,n_r]);
    
    y_ch = y_e - ( c_e .* h_e + b_e );
    
    nllh = nllh + 0.5*sum(sum(nansum(...
        ~isnan(y_e).*log(2*pi*noise_e) + y_ch.^2 ./noise_e,...
        3),2),1);
    
    if nargout > 1
        sh_e = permute(repmat(sim(ie).sy, [1,1,1,n_r]), [1,2,4,3]);
        dy_ch = - c_e .* sh_e;
        
        grad = grad + squeeze(sum(sum(nansum(...
            y_ch./noise_e .* dy_ch,...
            1),2),3));
        
        if nargout > 2
            for j=1:n_theta
                for k=1:n_theta
                    fim(j,k) = fim(j,k) + sum(sum(nansum(...
                        dy_ch(:,:,:,j) .* dy_ch(:,:,:,k) ./ noise_e)));
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

end


function [ varargout ] = nllh_laplace( sim, D, b, c, noise )

n_e = size(D,2);
if nargout > 1
    n_theta = size(sim(1).sy,3);
end

% initialize output
nllh = 0;
if nargout > 1
    grad = zeros(n_theta,1);
end

% compute output
for ie = 1:n_e
    b_e = b{ie};
    c_e = c{ie};
    noise_e = noise{ie};
    
    y_ch = bsxfun(@minus,D(ie).Y,bsxfun(@times,c_e,sim(ie).y)+b_e);
    
end

varargout{1} = nllh;
if nargout > 1
    varargout{2} = grad;
end

end

