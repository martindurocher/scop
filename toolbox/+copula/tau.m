function tau = tau(x)
%
% SUMMARY: 
%
%  Wrapper for the kendall tau function. See corr(x,'type','kendall')
%
%% CALL
%
%    tau = kendall_tau(x)
%
%  x: A  matrix with observation in rows (nx2).
%
%% AUTHOR: <martin.durocher@uqtr.ca>
    if(size(x,2) == 2)
        tmp = corr(x,'type','kendall');
        tau = tmp(2);
    else
        error('Must be bivariate observation');
    end

end
