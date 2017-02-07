function u = pbt2_copula(x,param)
%% WARNING: NOT TESTED
%
%% Summary
%
% Calculate the cdf of the bivariate t-squared copula 
%
%% CALL
%
% u = pbt_copula(w,param)
% 
% w : Matrix (nx2) of bivariate uniform variable
% param : [rho,nu,a] where rho = Correlation coefficient
%                          nu = degree of freedom
%                          a  =  non centrality parameter
%
%% Author <martin.durocher@uqtr.ca>
    
    if(param(3) == 0)
        u = cdf0(x,param);    
    else
        %u = cdf1(x,param);
        error('not implemented')
    end  
    
end

function u = cdf0(w,param)
    xp1 = (1+w(:,1))./2;
    xp2 = (1+w(:,2))./2;
    xm2 = (1-w(:,2))./2;
    u = 2* ( pbt_copula([xp1,xp2],param(1:2)) ...
             - pbt_copula([xp1,xm2],param(1:2)) ) - w(:,2);
end

function u = cdf1(w,param)
    w;param;
end