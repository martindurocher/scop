function z = rmvt(n,sigma,df)
%%
%% DESCRITION:
%
%    Simulating a multivariate Student vector using 
%       Singular Values Decomposition (SVD)
%% CALL
%
%    y = rmvnorm(n,sigma,mu)
%
%  n     : Size of the sample (Interger)
%  sigma : Covariance matrix (Matrix n x n )
%  mu    : Mean vector (default = 0)
%
%% NOTE
%
% see rmvnorm

if(df <= 0)
    error('Wrong')
end

w = sqrt(df./chi2rnd(df,1,n));
z = rmvnorm(n,sigma);

for(kk = 1:n)
    z(kk,:) = w(kk)*z(kk,:); 
end