function U = sim_multi(n,model,Sigma,param)
%%
%% SUMMARY:
%
%    Simulating a multivariate copula
%
%% CALL
%
%    U = sim_copula_multi(n,model,Sigma,param)
%
%  n     : Size of the sample (Interger)
%  model : Copula name (String)
%  Sigma : Correlation matrix (Matrix n x n )
%  param : Copula specific parameter (Numeric)
%
%% NOTE
%    
%  Avalaible copulas
%    Normal      : model = 'normal',  param = empty.
%    Student     : model = 'student', param = degree of freedom
%    Chi-squared : model = 'chisq',   param = non-centrality
%    t2-copula   : model = 't2',      param = [nu,a]
%    
%  In the case of the the smith model(model = 'smith'): 
%  The model is a spatial case of the Husler-Reiss copula .
%  Sigma = coordinate of the simulted point,
%  param = a 2x2 covariance matrix.
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>
%  adaptation of original codes from JF Quessy.

    d = size(Sigma,1);
    mu = zeros(1,d);
    
    switch model
        case 'normal'
            X = rmvnorm(n,Sigma);
            U = normcdf(X);
        
        case {'student','t'}
            X = mvtrnd(Sigma,param(1),n);
            U = tcdf(X,param(1));
        
        case 'chisq'
            Y = abs(rmvnorm(n,Sigma)+param(1));
            if(param(1) == 0)
                U = 2*normcdf(Y) - 1;
            else
                U = normcdf(Y-param(1)) + normcdf(Y+param(1)) - 1;
            end
            
        case 'smith'
            U = gevcdf(rsmith(n,Sigma,param),1,1,1);
            
        case 't2'
            Y = abs(mvtrnd(Sigma,param(1),n)+param(2));
            if(param(2) == 0)
                U = 2 * tcdf(Y,param(1)) - 1;
            else
                U = tcdf(Y-param(2),param(1)) + ...
                    tcdf(Y+param(2),param(1)) - 1;
            end
            
        otherwise
            error('Invalid choice of copula');
        
    end
end
    

