function pdf = dbchisq(u,rho,a,uselog, tab)
%%
%% SUMMARY:
%
%    Calculate the Density of a bivariate chi-squared copula
%
%% CALL 
%
%    pdf = dbchisq_copula(u,rho,a,uselog)
%  
%  pdf    : return density (Numeric)
%  u      : Uniform sample (Matrix n x 2)
%  rho    : Correlation coefficient (Numeric)
%  a      : Copula specific parameter (Numeric)
%  uselog : Should the log density be return (logic)
%  tab    : Tabulated value of the ha function (if necessary)
% 
%% Author: <martin.durocher@uqtr.ca>

    TOL_A = 1e-6;

    if(nargin < 4); uselog = 0;end
    if(nargin < 5); tab = [];end
        
    if( a < TOL_A)
        up = (1 + u(:,1))/2; 
        un = (1 - u(:,1))/2; 
        v  = (1 + u(:,2))/2;
    
        pdf = 0.5 *(copula.dbnorm([up,v],rho) + ...
                                    copula.dbnorm([un,v],rho));
    
    else
        % if a table of value is provided for interpolation
        if(~isempty(tab))
            ha = interp1(tab(:,1),tab(:,2),u,'spline');
            dha = interp1(tab(:,1),tab(:,3),u,'spline');
        else
           [ha,dha] = copula.chisq_ha(u,a);
        end
        
        uha = u + ha;
        mdha = 1 - dha;
        
        pdf = mdha(:,1) .* mdha(:,2) .* copula.dbnorm(uha,rho) + ...
              mdha(:,1) .*  dha(:,2) .* copula.dbnorm([uha(:,1),ha(:,2)],rho) + ...
               dha(:,1) .* mdha(:,2) .* copula.dbnorm([ha(:,1),uha(:,2)],rho) + ...
               dha(:,1) .*  dha(:,2) .* copula.dbnorm(ha,rho);
        
    end 
    
    if(uselog)
        pdf=log(pdf);
    end
end
    
    
