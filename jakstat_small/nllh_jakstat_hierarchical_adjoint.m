function [ varargout ] = nllh_jakstat_hierarchical_adjoint(theta,D,scOptions)

amiOptions.rtol = 1e-12;
amiOptions.atol = 1e-14;
simfun = @simulate_jakstat_hierarchical_adjoint;

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
% ne = size(D,2);
% 
% % forward simulation
% 
% sim = struct([]);
% for ie = 1:ne
%     
%     kappa_e = [kappa(:,ie); 1; 1];
%     
%     sol = simulate_jakstat_hierarchical_adjoint(D(ie).t,theta,kappa_e,[],amiOptions);
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
% [b,c,sigma2,b_by_y,c_by_y,sigma2_by_y] = opt_scalings(sim,D,scOptions);
% 
% % nllh,snllh, s2nllh
% 
% % prepare sensi
% switch nargout
%     case 1
%         amiOptions.sensi = 0;
%     case 2
%         amiOptions.sensi = 1;
%     otherwise
%         amiOptions.sensi = 2;
% end
% 
% if nargout == 1
%     nllh = hieropt_nllh_forward(false,sim,D,b,c,sigma2);
% else
%     nllh = 0;
%     snllh = 0;
%     if nargout > 2
%         s2nllh = 0;
%     end
%     
%     for ie = 1:ne
%         nr = size(D(ie).my,3);
%         for ir = 1:nr
%             clear amiData
%             
%             c_re = c_by_y{ie}(1,:,ir);
%             c_re = reshape(c_re,[],1);
%             
%             sigma2_re = sigma2_by_y{ie}(1,:,ir);
%             sigma2_re = reshape(sigma2_re,1,[]);
%             sigma2_re = repmat(sigma2_re,length(D(ie).t),1);
%             
%             kappa_re = [kappa(:,ie); c_re(1:2)];
%             
%             amiData.t = D(ie).t;
%             amiData.Y = D(ie).my(:,:,ir);
%             amiData.Sigma_Y = sqrt(sigma2_re);
%             amiData.condition = kappa_re;
%             amiData = amidata(amiData);
%             
%             sol = simulate_jakstat_hierarchical_adjoint(D(ie).t,theta,kappa_re,amiData,amiOptions);
%             
%             if sol.status ~= 0
%                 error('Could not integrate ODE.');
%             end
%             
%             nllh = nllh - sol.llh;
%             snllh = snllh - sol.sllh;
%             if nargout > 2
%                 s2nllh = s2nllh - sol.s2llh;
%             end
%         end
%     end
% %     fprintf('%.15f %.15f %.15f\n',varargout{1},cost,llh);
% %     disp(mat2str(varargout{2}));
% end
% 
% varargout{1} = nllh;
% if nargout > 1
%     varargout{2} = snllh;
%     if nargout > 2
%         varargout{3} = s2nllh;
%     end
% end
% 
% end