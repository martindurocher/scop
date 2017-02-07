function [ha,dha]= chisq_ha(u,a)
%%
%% SUMMARY:
%
%    Compute an intermediate value use in the Density of a 
%     bivariate chi-squared copula. Can be passed as a table 
%     in spmodel to speed up the estimation of a chi-squared 
%     copula.
%
%% CALL 
%
%    [ha, dha] = copula.chisq_ha(u,a)
%  
%  u : Matrix of pseudo-observations in [0,1];
%  a : Copula specific parameter of the Chi-squared copula (Scalar).
%
%% Author: <martin.durocher@uqtr.ca>

    ipsi = psi_inv(u ,a) + a;
    ha = 1-normcdf(ipsi);

    if(nargout > 1)
    dha = normpdf(ipsi) ./(normpdf(ipsi) + normpdf(ipsi - 2*a));
    end

end

function y = psi_inv(u,a) 
    y = zeros(size(u));
    
    for i=1:size(u,1)
        for j=1:size(u,2)
            y(i,j) = fzero(@(x)(normcdf(x - a) + normcdf(x + a) - 1 ...
                - u(i,j)),[0,11]);
        end
    end
    
end