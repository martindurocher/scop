function u = rbhusler(n,a, method)
%
%% SUMMARY
%
% Generate bivariate sample from the Husler-Reiss Copula.
%
%% CALL
%
%    u = rbhusler_copula(n,a, method)
%   
% n      : number of simulated value
% a      : Dependence parameter (a>0)
% method : Simulation method
%
%% DETAILS
%
% The method can be obtained either by the conditional copula ('cond')
% or approximated by using two stations of a Smith model ('smith')  

%% <martin.durocher@uqtr.ca>

if(nargin <3); method = 'cond'; end

switch method
    case 'cond'
        u = sim_cond(n,a);
        
    case 'smith'
        u = sim_smith(n,a);
   
end % switch
end % function

 function u = sim_cond(n,a)
    u = rand(n,2);
    POSMIN = 1e-6;
    for(j = 1:n)
        u(j,2) = fzero(@(v) ifun(v,u(j,2),u(j,1),a),[POSMIN,1-POSMIN]); 
    end
 end

 function u = sim_smith(n,a)
    sp = spmodel('smith','iso',a);
    u = sp.sim_coord(n,[0,0;0,1]);
 end
 

function d = ifun(v,p,u,a)
    x = -log(u); y = -log(v);
    lxy = log(x) - log(y);
    ia = 1./a;
    a2 = a./2;
    nx = normcdf(ia + a2 .* lxy);
    ny = normcdf(ia - a2 .* lxy);
    d = x .* (1 - nx) - y .* ny + log(nx) - log(p);
end
    