function u = pbnorm(w,rho)
%
%% Summary
%
% Calculate the bivariate normal copula
%
%% CALL
% u = pbnorm_copula(w,rho)
% 
% w : Matrix (nx2) of bivariate uniform variable
% rho : Correlation coefficient
%
%% Author <martin.durocher@uqtr.ca>

if(size(w,2) ~= 2) 
     error('Not bivariate vectors')
end
  
y = norminv(w);
u = mvncdf(y,[0,0],[1,rho;rho,1]);

end