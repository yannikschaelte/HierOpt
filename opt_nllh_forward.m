function [ varargout ] = opt_nllh_forward( sim, D, b, c, sigma2 )
% compute the negative log-likelihood (sum of weighted squares), typically
% for minimization
% 
% Input:
%
% Output:
% * nllh
% * grad
% * fim
%
% History:
% 2018/01/10 Yannik Schaelte

ne = size(D,2);
if nargout > 1
    n_theta = size(sim(1).sy,3);
end

nllh = 0;
if nargout > 1
    grad = zeros(n_theta,1);
    if nargout > 2
        fim = zeros(n_theta);
    end
end


for ie = 1:ne
    nr = size(D(ie).my,3);
    
    b_e = b{ie};
    c_e = c{ie};
    sigma2_e = sigma2{ie};
    
    y_ch = bsxfun(@minus,D(ie).my,bsxfun(@times,c_e,sim(ie).y)+b_e);
    
    nllh = nllh + 0.5*sum(sum(nansum(bsxfun(@times,~isnan(D(ie).my),log(2*pi*sigma2_e))+...
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
                    fim(j,k) = fim(j,k) + nansum(nansum(nansum(bsxfun(@rdivide,...
                        bsxfun(@times,dy_ch(:,:,:,j),dy_ch(:,:,:,k)),sigma2_e))));
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

