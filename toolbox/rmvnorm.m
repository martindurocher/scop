function y=rmvnorm(n,sigma,mu)
%%
%% DESCRITION:
%
%    Simulating a multivariate Normal vector using 
%       Singular Values Decomposition (SVD)
%% CALL
%
%    y = rmvnorm(n,sigma,mu)
%
%  n     : Size of the sample (Interger)
%  sigma : Covariance matrix (Matrix n x n )
%  mu    : Mean vector (default = 0)
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>

% verify that sigma is symetric 
if(min(min(sigma ~= sigma')))
    error('The matrix is not symetric')   
end

% extract the vector dimension
m = size(sigma,1);

% set default mean vector if necessary
if(nargin < 3) mu = zeros(1,m); end

% simulate a independant sample
y = randn(n,m);

%Decompose the covariance matrix by SVD
[u,s,v] = svd(sigma);

% return a transformed vector
y = repmat(mu,n,1)+ y * sqrt(s) * v' ;