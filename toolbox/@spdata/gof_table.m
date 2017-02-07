function [out,sp0] = gof_table(obj, type, nboot, estim)

% to finish
if(nargin < 2); type  = 'all'; end
if(nargin < 3); nboot = 1000; end
if(nargin < 4); estim = 'pairwise'; end
    
[out,sp0] = obj.gof(nboot,type,estim,'table');

 