function obj = set_ha(obj,n)
%% Method
%
% Create a table of intermediate values that are used to
% speed up the evaluation of the likelihood of a chi-squared
% copula. See copula.chisq_ha
%
%% CALL
%
%  obj = set_ha(obj,n)
%
%  n = number of point for interpolation (default = 500)
%
%% NOTE
%
% Use n=0 to remove the table.
%%%%%%%%%%%%%%%%%%%%%%%%%

if(~strcmp(obj.mdl,'chisq'))
    error('Must be a Chi-squared copula')
end

%default value
if(nargin <3); n = 500; end

if(n == 0)
    obj.tab_ha = [];
else
    tab_u = linspace(0,1,n);
    [tab_ha,tab_dha] = copula.chisq_ha(tab_u,obj.a);
    obj.tab_ha = [tab_u;tab_ha;tab_dha]';
end