function theta = tau_theta(tau,model,param)
%%
%% SUMMARY:
%
%    Returns the parameter of a copula from Kendall's
%
%% CALL
%
%    theta = tau_theta(tau,model,param)
%
%  theta : Correlation parameter of the copula
%  tau   : Kendall's tau (Matrix)
%  model : Copula name (String)
%  param : Copula specific parameter (Scalar)
%
%% NOTE
%    
%  Avalaible copulas
%    Normal      : model = 'normal',  param = empty.
%    Student     : model = 'student', param = degree of freedom
%    Chi-squared : model = 'chisq',   param = non-centrality
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>
error('Not yet')
    % Numerical approximation of other cases
    MAXPOSCORR = 0.999999 ;
    MINPOSCORR = 0.000001 ;
    
    if(strcmp(model,'normal'))
        theta = sin((pi/2)*tau);
        
    elseif(strcmp(model,'student'))
        theta = sin((pi/2)*tau);
        
    elseif(strcmp(model,'chisq'))
        [m,n] =  size(tau);
        theta = zeros(m,n);
        
        for i = 1:m
            for j = 1:n
                if(tau(i,j) == 0)
                    theta(i,j) = 0;
                elseif(tau(i,j) == 1)
                    theta(i,j) = 1;
                elseif(tau(i,j) == -1)
                    theta(i,j) = -1;
                else
                    fun = @(x) fun_tau_chisq(x,tau(i,j),param);
        
                    if(tau < 0)
                        theta(i,j) = fminbnd(fun,-MAXPOSCORR,MINPOSCORR);
                    else
                        theta(i,j) = fminbnd(fun,MINPOSCORR,MAXPOSCORR);
                    end
                end
                
            end
        end
        
    else
        error('Invalid family of copula');
    end
end
	
% Intermediate function to solve 
function out = fun_tau_chisq(c,tau0,a)
    s2a = sqrt(2)*a;
    mu = [0,0];
    sig = [1,c;c,1];
    
    z = (4*mvncdf([s2a,s2a],mu,sig)- 4*normcdf(s2a) + 1);
    tau_chisq = 2 / pi * asin(c) * z;
    out = (tau0 - tau_chisq)^2;

end
 