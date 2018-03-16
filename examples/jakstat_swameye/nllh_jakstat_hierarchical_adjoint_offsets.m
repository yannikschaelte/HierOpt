function [ varargout ] = nllh_jakstat_hierarchical_adjoint_offsets(theta,D,scOptions)

amiOptions.rtol = 1e-12;
amiOptions.atol = 1e-14;
simfun = @simulate_jakstat_hierarchical_adjoint_offsets;

switch nargout
    case 1
        [varargout{1}] = hieropt_nllh_adjoint(simfun,theta,D,amiOptions,scOptions);
    case 2
        [varargout{1},varargout{2}] = hieropt_nllh_adjoint(simfun,theta,D,amiOptions,scOptions);
    case 3
        [varargout{1},varargout{2},varargout{3}] = hieropt_nllh_adjoint(simfun,theta,D,amiOptions,scOptions);
end

end

% amiOptions.sensi_meth = 'adjoint';
% amiOptions.sensi = 0;
% 
% % for every experiment and replicate, do the optimization
% n_e = size(D,2);
% 
% % forward simulation
% 
% sim = struct([]);
% for ie = 1:n_e
%     
%     kappa_e = [kappa(:,ie); 0; 0; 1; 1];
%     
%     sol = simulate_jakstat_hierarchical_adjoint_offsets(D(ie).t,theta,kappa_e,[],amiOptions);
%     
%     if (sol.status ~= 0)
%         error('Could not integrate ODE.');
%     end
%     
%     sim(ie).y = sol.y;
% end
% 
% % optimal scalings
% 
% [c,sigma2,c_by_y,sigma2_by_y,b,b_by_y] = opt_scalings_normal(sim,D,scOptions);
% 
% % nllh,snllh, s2nllh
% 
% % prepare output and sensi
% switch nargout
%     case 1
%         varargout{1} = 0;
%         amiOptions.sensi = 0;
%     case 2
%         varargout{1} = 0;
%         varargout{2} = 0;
%         amiOptions.sensi = 1;
%     case 3
%         varargout{1} = 0;
%         varargout{2} = 0;
%         varargout{3} = 0;
%         amiOptions.sensi = 2;
%     otherwise
%         error('Only supports up to 3 outputs.');
% end
% 
% if nargout == 1
% llh = 0;
%     for ie = 1:n_e
%         sigma2_e = sigma2(:,:,:,ie);
%         c_e = c(:,:,:,ie);
%         b_e = b(:,:,:,ie);
%         y_ch = bsxfun(@minus,D(ie).my,bsxfun(@add,bsxfun(@times,c_e,sim(ie).y),b_e));
%         llh = llh + 0.5*sum(sum(nansum(bsxfun(@times,~isnan(D(ie).my),log(2*pi*sigma2_e))+...
%             bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2_e),1),3),2);
%     end
% else
%     for ie = 1:n_e
%         n_r = size(D(ie).my,3);
%         for ir = 1:n_r
%             clear amiData
%             
%             c_re = c_by_y(:,:,ir,ie);
%             c_re = reshape(c_re,[],1);
%             
%             b_re = b_by_y(:,:,ir,ie);
%             b_re = reshape(b_re,[],1);
%             
%             sigma2_re = sigma2_by_y(:,:,ir,ie);
%             sigma2_re = reshape(sigma2_re,1,[]);
%             sigma2_re = repmat(sigma2_re,length(D(ie).t),1);
%             
%             kappa_re = [kappa(:,ie); b_re(1:2); c_re(1:2)];
%             
%             amiData.t = D(ie).t;
%             amiData.Y = D(ie).my(:,:,ir);
%             amiData.Sigma_Y = sqrt(sigma2_re);
%             amiData.condition = kappa_re;
%             amiData = amidata(amiData);
%             
%             sol = simulate_jakstat_hierarchical_adjoint_offsets(D(ie).t,theta,kappa_re,amiData,amiOptions);
%             
%             if sol.status ~= 0
%                 error('Could not integrate ODE.');
%             end
%             
%             varargout{1} = varargout{1} - sol.llh;
%             varargout{2} = varargout{2} - sol.sllh;
%             if nargout > 2
%                 varargout{3} = varargout{3} - sol.s2llh;
%             end
%         end
%     end
% %     fprintf('%.15f %.15f %.15f\n',varargout{1},cost,llh);
% %     disp(mat2str(varargout{2}));
% end
% 
% end