function pdf = dbt2(u,rho,nu,a,uselog,isu)
%% WARNING NOT TESTED
%
%% SUMMARY:
%
%    Calculate the Density of bivariate t-squared copula
%
%% CALL
%
%  pdf = dbt2_copula(u,param,uselog, tol)
%  
%  u      : Uniform sample (Matrix n x 2)
%  param  : vector of [rho,nu,a] respectively the
%           the coefficient of correlation, the degree of freedom
%           and the non-centrality parameter
%  uselog : return log density {0,1} (default = 0, FALSE)
%  tol    : tolerance value to use the specific equation for a=0
%
%% Author: <martin.durocher@uqtr.ca>

if(nargin < 5); uselog = 0; end
if(nargin < 6); isu = 1; end % flad to avoid repetition of tinv
        
if( a < 1e-8)
    if(isu)
        up = (1 + u(:,1))/2; 
        un = (1 - u(:,1))/2; 
        v  = (1 + u(:,2))/2;
    
        pdf = 0.5 * (copula.dbt([up,v],rho,nu) + ...
                     copula.dbt([un,v],rho,nu));
    else
        up = u(:,1); 
        un = u(:,2); 
        v  = u(:,3);
    
        pdf = 0.5 * (copula.dbt([up,v],rho,nu,0,0) + ...
                     copula.dbt([un,v],rho,nu,0,0));
    end
        
    
else
    error('not implemeted yet')
    
end

if(uselog)
    pdf = log(pdf);
end
