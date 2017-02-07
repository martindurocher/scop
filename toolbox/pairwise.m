function y = pairwise(x,opt)
%
%% DESCRIPTION
%
%  Create a matrix of all the bivariate pairs.
%
%% CALL
%
%   y = pairwise(x)
%
% x: Row matrix (n x 1)
% y: Row matrix [m(m-1)/2 x 2n] of the pairs of rows 
%
%
%% Author: martin du rocher <martin.durocher@uqtr.ca>

if(nargin < 2); opt ='d'; end; 

[m,n] = size(x);
mpair = m*(m-1)/2;

if(m==1)
    y=x;
elseif(n==1)
    y = src.pairwise_c(x);
else
    y = zeros(mpair,2*n);
    for i=1:n
        z = src.pairwise_c(x(:,i));
        y(:,i) = z(:,1);
        y(:,i+n)= z(:,2);
    end
    
end

switch opt
    case {'p','permute','perm'}
        for(j=1:mpair)
            if(rand(1)>.5)
                y(j,:) = [y(j,2),y(j,1)];
            end
        end
end