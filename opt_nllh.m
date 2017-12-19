function [ nllh, snllh ] = opt_nllh( D, sim, b, c, sigma2 )
% compute the negative log-likelihood (sum of weighted squares), typically
% for minimization

n_e = size(D,2);

nllh = 0;
if nargout > 1
    snllh = 0;
end

for ie = 1:n_e
    n_r = size(D(ie).my,3);
    sigma2_e = sigma2(:,:,:,ie);
    c_e = c(:,:,:,ie);
    b_e = b(:,:,:,ie);
    y_ch = bsxfun(@minus,D(ie).my,bsxfun(@times,c_e,sim(ie).y)+b_e);
    nllh = nllh + 0.5*sum(sum(nansum(bsxfun(@times,~isnan(D(ie).my),log(2*pi*sigma2_e))+...
        bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2_e),1),3),2);
    if nargout > 1
        dy_ch = - bsxfun(@times,permute(c_e,[1,2,4,3]),repmat(sim(ie).sy,[1,1,1,n_r]));
        snllh = snllh + permute(sum(sum(nansum(...
            bsxfun(@times,permute(...
            bsxfun(@rdivide,y_ch,sigma2_e),[1,2,4,3]),...
            dy_ch),1),4),2),[3,2,1]);
    end
end

end

