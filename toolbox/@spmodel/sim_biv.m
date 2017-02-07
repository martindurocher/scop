function out = sim_biv(obj,nobs,rho, method)
%% Method: sim_biv
%
% SUMMARY: Simulate a bivariate copula. 
%
% CALL: out = spmodel.sim_biv(nobs, rho, method)
%
%   nobs  : number of simulations
%   rho   : dependance coefficient of the bivariate copula
%   method: method used in case of a smith model 
%           see rsmith.m
%
% NOTE: see also sim_copula_multi.m and rsmith.m
%   

switch obj.mdl
    case {'normal','chisq','student','t2'}
        sig = [1,rho;rho,1];
        tmp = obj.theta();
        out = copula.sim_multi(nobs, obj.mdl,sig,tmp(1:2));
    case 'smith'
        if(nargin < 4); method = 'smith';end
        out = copula.rbhusler(nobs,rho,method);
    case 'ind'
        out = rand(nobs,2);
end