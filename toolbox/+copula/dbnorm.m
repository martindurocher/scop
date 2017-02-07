function pdf = dbnorm(u,rho, uselog)
%%
%% SUMMARY:
%
%    Calculate the Density of bivariate normal copula
%
%% CALL
%
%    pdf = dbnorm_copula(u,rho, uselog)
%  
%  pdf   : returned density (Numeric)
%  u     : Uniform sample (Matrix n x 2)
%  rho   : Correlation coefficient (Vector)
%  uselog : If the log-denstiy should be return
%
%% Author: <martin.durocher@uqtr.ca>

    if(nargin<3); uselog = 0; end
    
    if(size(rho,1) == 1 && size(rho,2)> 1)
        rho = rho';
    end
    
    x = norminv(u(:,1),0,1);
    y = norminv(u(:,2),0,1);
    r2 = 1 - (rho .* rho);
    pdf = log(1/2/pi) - 0.5 * log(r2) ...
          - ((x.* x) - 2 * (x .* y .* rho) + (y .* y)) ./ (2 * r2)...
          - log(normpdf(x,0,1)) - log(normpdf(y,0,1));
    
    if(not(uselog))
        pdf=exp(pdf);
    end
end