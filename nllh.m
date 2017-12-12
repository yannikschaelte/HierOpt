function [ cost, gradient ] = nllh( data, sim, b, c, sigma2 )
% compute the negative log-likelihood (sum of weighted squares), typically
% for minimization
cost = 0;
for ie = 1:n_e
    sigma2_e = sigma2(:,:,:,ie);
    c_e = c(:,:,:,ie);
    y_ch = bsxfun(@minus,D(ie).my,bsxfun(@times,c_e,sim(ie).y));
    nlLH = nlLH + sum(sum(nansum(bsxfun(@times,~isnan(D(ie).my),log(2*pi*sigma2_e))+...
        bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2_e),1),3),2);
end

end

