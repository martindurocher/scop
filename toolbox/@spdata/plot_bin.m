function plot_bin(obj,type)
% SUMMARY : Plot a coefficient of the group of lag distances
%
% CALL: obj = obj.plot_bin(type)
%
%    type: 'rho' = Spearman rho (default)
%          'tau' = kendall tau

if(nargin < 2) type = 'rho'; end
         
switch type
    case 'rho'
        scatter([0,obj.mids],[1,obj.bin_rho()]);
    case 'tau'
        scatter([0,obj.mids],[1,obj.bin_tau()]);
end
