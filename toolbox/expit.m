function y = expit(x)
%%
%% DESCRIPTION:
%
%    Calculate the logit function.
%
%% CALL
%
%    y = expit(x)
%  
%  x : Range of the link function (Numeric)
%
%% NOTE
%
% see also expit its inverse
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>

id1 = (x==Inf);
id2 = (x~=Inf);
y(id1) = 1;
y(id2) = exp(x(id2))./(1+exp(x(id2)));
