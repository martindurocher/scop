function V = transform(X,method)  
%% SUMMARY 
%
%  Compute the probability integral transform of a bivariate sample
%
% CALL
%
% V = kendall_transform(X,method)
%  
% x      : sample (Matrix n x 2)
% method : method to be used ('char')
%
%% NOTE
%
% There is two method: One that use excheangeability ('xch') and 
%  the one without ('normal'). By defaut 'normal'.
%
% Author: Martin du Rocher <martin.durocher@uqtr.ca>
    
if(nargin < 2); method = 'normal'; end
if(~all_isfinite(X)); error('Non finite values'); end
if(size(X,2) ~= 2);   error('Wrong dimension'); end
    
switch method
    case {'normal','cvm'}
        V = src.kendall_transform_c(X); 
    case 'xch'
        V = src.kendallx_transform_c(X);
end

