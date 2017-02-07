function y = logit(x)
%
%% DESCRIPTION:
%
%    Calculate the logit function.
%
%% CALL
%
%    y = logit(x)
%  
%  x : Range of the link function (Numeric)
%
%% NOTE
%
% see also expit its inverse
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>

y = log(x ./(1-x));