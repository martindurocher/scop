function L = get_lik(obj,param)
%% Method: get_lik
%
% SUMMARY : evaluate the likelihood for the spmodel
%
% CALL: L = spdata.get_lik(param)
%
%   param : parameter to evaluate likelihood  
%    
switch obj.sp.mdl
    case 'normal'
    
        zdata = norminv(obj.ranks)';
        mu = zeros(1,size(zdata,1));
        
        nparam = size(param,1);
        for(ii = 1:nparam)
            sigma = cor_model(obj.distance, ...
                obj.sp.link,param(ii,:));
                                
            L(ii) = mvnpdf(zdata,mu,sigma);
        end
        
     case 'student'
        zdata = tinv(obj.ranks,obj.sp.df)';
        mu = zeros(1,size(zdata,1));
        
        nparam = size(param,1);
        for(ii = 1:nparam)
            sigma = cor_model(obj.distance, ...
                obj.sp.link,param(ii,:));
            
            L(ii) = mvnpdf(zdata,mu,sigma);
        end
            
    otherwise
        error('Not implememented yet')
end
