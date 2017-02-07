function u = pbt(w,param)
%
%% Summary
%
% Calculate the cdf of the bivariate t-copula 
%
%% CALL
% u = pbt_copula(w,param)
% 
% w : Matrix (nx2) of bivariate uniform variable
% param : [rho,nu] where rho = Correlation coefficient
%                         nu = degree of freedom
%
%% Author <martin.durocher@uqtr.ca>
if(size(w,2) ~= 2) 
     error('Not bivariate vectors')
end

y = tinv(w,param(2));
u = mvtcdf(y,[1,param(1);param(1),1],param(2));
end