function V = all_isfinite(x)
%% SUMMARY:
%
%    Test if all element of a matix are finite. See isfinite.m
%
%% CALL
%
%    V = all_isfinite(x)
%
%  x : A matrix
%  V : A logical {0,1} 
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>
V = min(min(isfinite(x)));