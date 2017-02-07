function tau = theta_tau(theta,model,param)
%%
%% SUMMARY:
%
%    Calculate the Kendall's tau base on the copula parameter
%
%% CALL
%
%   tau = theta_tau(theta,model,param)
%
%  tau   : Kendall's tau (Matrix)
%  theta : Correlation parameter of the copula
%  model : Copula name (String)
%  param : Copula specific parameter (Scalar)
%
%% NOTE
%    
%  Avalaible copulas
%    Normal      : model = 'normal',  param = empty.
%    Student     : model = 'student', param = degree of freedom
%    Chi-squared : model = 'chisq',   param = non-centrality

%% Author: Martin du Rocher <martin.durocher@uqtr.ca>
     
     
     [m,n] = size(theta);
     
    switch model
        case {'normal', 'student'}
            tau =  2 / pi * asin(theta); 
        case 't2'
            error('not implemented yet')
        case 'chisq'
            
            tau = zeros(m,n);
            mu = [0,0];
            s2a = sqrt(2)*param; 
            
            for i=1:m
                for j = 1:n
                    if(theta(i,j) == 1)
                        tau(i,j) = 1;
                    elseif(theta(i,j) == -1)
                        tau(i,j) = -1;
                    elseif(theta(i,j) == 0)
                        tau(i,j) = 0;
                    else
                        sig = [1,theta(i,j);theta(i,j),1];  
                        z = (4*mvncdf([s2a,s2a],mu,sig)- ...
                        4*normcdf(s2a) + 1);
                        tau(i,j) = 2 / pi * asin(theta(i,j)) * z; 
                    end
                end
            end
  
        otherwise
            error('Invalid family of copula');
    end
    
end

function dtau = smith_B(w,a)
    mw = 1-w;
    ia = 1/a;
    a2 = a/2 ;
    nw  = normcdf(ia + a2 .* log(w./mw));
    n1w = normcdf(ia + a2 .* log(mw./w));
    
    Bprime = (nw - n1w);
    B = w .* nw + mw .* n1w; 
    dtau = (2*w-1).* Bprime .* B + w .* mw .* Bprime .* Bprime;
    dtau = dtau ./ (B .* B);
end
 