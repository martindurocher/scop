function pdf = dbt(u,rho,nu,uselog,isu)
%%
%% SUMMARY:
%
%    Calculate the Density of bivariate t-copula
%
%% CALL
%
%  pdf = dbt_copula(u,param,uselog, tol)
%  
%  u      : Uniform sample (Matrix n x 2)
%  param  : vector of [rho,nu,a] respectively the
%           the coefficient of correlation, the degree of freedom
%           and the non-centrality parameter
%  uselog : return log density {0,1} (default = 0, FALSE)
%
%% NOTE

%default parameter
if(nargin<4); uselog = 0; end
if(nargin<5); isu = 1; end % flag to repetitive call of tinv
    
if(size(rho,1) == 1 && size(rho,2)> 1)
    rho = rho';
end

if(isu)
    x = tinv(u(:,1),nu);
    y = tinv(u(:,2),nu);
else
    x = u(:,1);
    y = u(:,2);
end

detr = 1-rho.^2;

cst = gamma(nu/2+1) ./ (gamma(nu/2) * pi * nu * sqrt(detr));
pdf = 1 + (x.^2 - 2 * x .* y .* rho + y.^2) ./ (detr * nu);
pdf = cst .* (pdf).^(-(nu/2+1)) ./ tpdf(x,nu) ./ tpdf(y,nu);

if(uselog)
    pdf = log(pdf);
end