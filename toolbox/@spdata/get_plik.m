function L = get_plik(obj,param)
%% Method: get_plik
%
% SUMMARY : evaluate the pairwise likelihood
%
% CALL: L = spdata.get_plik(param)
%
%   param : set of parameters at which to evaluate the likelihood
%   h     : Pairwise distances (matrix 1 x n)
%   u     : Pairs of sites (matrix n x 2)   
%

h = obj.pdistance;
u = obj.pairs;

switch obj.sp.mdl
    case 'normal'
        nparam = size(param,1);
        for(ii = 1:nparam)
            rho = cor_model(h,obj.sp.link,param(ii,:));
            L(ii) = sum(dbnorm_copula(u,rho,1));
        end
                
    case 'chisq'
        nparam = size(param,1);
        for(ii = 1:nparam)
            rho = cor_model(h,obj.sp.link,param(ii,2:4));
            L(ii) = sum(dbchisq_copula(u,rho, param(ii,1),1));
        end
        
    case 'student'
        nparam = size(param,1);
        for(ii = 1:nparam)
            rho = cor_model(h,obj.sp.link,param(ii,2:4));
            L(ii) = sum(dbt_copula(u,rho,param(ii,1),1));
        end
        
    case 't2'
        nparam = size(param,1);
        for(ii = 1:nparam)
            rho = cor_model(h,obj.sp.link,param(ii,3:5));
            L(ii) = sum(dbt2_copula(u,rho, ...
                            param(ii,1),param(ii,2),1));
        end
        
    otherwise
        error('Not implememented yet')
end

L = exp(L);
