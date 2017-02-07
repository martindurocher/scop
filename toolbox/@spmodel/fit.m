function obj = fit(obj,method,h,u)
%% Method: fit
%
% SUMMARY: Estimate the parameter of the spmodel. 
%
% CALL: spmodel = spmodel.fit(type,h,u)
%  
%   method: Estimation method (now only 'pairwise') 
%   h     : Pairwise distances (matrix 1 x n)
%   u     : Pairs of sites (matrix n x 2)
%
% NOTE: see also copula.fit_pairwise.m
%

if(strcmp(obj.mdl,'ind')) 
    error('Independent copula does not have parameters')
end

if(nargin < 2); method='pairwise';end
    
obj.nllik_type = method;

switch method
    case 'full'
        [theta, nllik] = copula.fit_full(h,u,obj);
    case 'pairwise'
        [theta, nllik] = copula.fit_pairwise(h,u,obj);
        %obj.nllik = nllik;
        %obj = obj.update(theta);
    otherwise
            error('Only pairwise is implemented');
end
   
obj.nllik = nllik;
obj = obj.update(theta);