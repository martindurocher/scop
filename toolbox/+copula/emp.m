function [ecdf,b1,b2] = emp(x,method,b1,b2)
%
%% SUMMARY
%
% Calculate the empirical copula of a bivariate sample. 
%   
%% CALL
%
%  [ecdf,b1,b2] = copula.emp(x,method,b1,b2)
%
% ecdf    : Distance criterion (return)
% b1,b2   : Breaks points for the respective variates (Matrix n x 1)
% x       : sample (Matrix n x 2)
% method  : Method to calculate the ecdf (char)
%
%% NOTE
%
%  There is two methods: one using exchangeability ('xch') and 
%   the other without ('normal'). By defaut 'normal' is used.
%
%  If no breaks are provided, the sampled values are used.
%
%
%% EXEMPLE 
%
%  u = rand(100,2);
%  bb1 = (1:10)/11;
%  bb2 = (1:20)/21;
%  [a,b,c] = emp_copula(u);
%  [a,b,c] = emp_copula(u,'xch',bb1,bb2);
%  
%% Author Martin du Rocher <martin.durocher@uqtr.ca>

if(nargin < 2); method = 'normal'; end
if(nargin < 3); b1 = unique(sort(x(:,1)')); end
if(nargin < 4); b2 = unique(sort(x(:,2)')); end

x = sortrows(x);

switch method
    case 'normal'
        ecdf = src.emp_copula_c(x,b1,b2);
    case 'xch'
        ecdf = src.empx_copula_c(x,b1,b2);
    otherwise
        error('Wrong method')
end