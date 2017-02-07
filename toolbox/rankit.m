function u = rankit(x)
 %
%% DESCRIPTION
%
%  Transform the data in pseudo observation, using n+1 as denominator (weibull)
%
%% CALL
%
%   u = rankit(x)
%
% x: Row matrix (m x n), input data
% u: Row matrix (m x n) output data in [0,1]
%
%
%% Author: martin du rocher <martin.durocher@uqtr.ca>   

[m,n] = size(x);
    
if(m == 1 && n >1)
    u = tiedrank(x')/ (n+1);
    u = u';
else
    u = tiedrank(x)/ (m+1);
end