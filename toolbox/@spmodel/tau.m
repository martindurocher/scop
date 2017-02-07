function out = tau(obj,x)
%% Method: tau
%
% SUMMARY: Return the kendall tau in respect of the distance. 
%  Except for the Smith model it correspond to a correlation coefficient
%
% CALL: out = spmodel.tau(x)
%           
% x : Distance where the 
%

switch obj.mdl
    case 'ind'
        out = zeros(size(x));
    case {'normal','student'}
        out = copula.theta_tau(rho(obj,x),obj.mdl);
    case {'chisq'}
        out = copula.theta_tau(rho(obj,x),obj.mdl,obj.a);
    otherwise
        error('not implemented')
end
