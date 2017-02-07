function cdf = pbhusler(u,a)
%
%% SUMMARY
%
% Cumulative distribution function of the Husler-Reiss Copula
%
%% CALL
%
%    cdf = pbhusler_copula(u,v,a)
%
% u : data (Matrix nx2)
% a : Dependence parameter (a>0) 
%
%% <martin.durocher@uqtr.ca>

v = u(:,2);
u = v(:,1);

% Transpose if necessary %
[m,n] = size(a)
if(m>1 & n == 1); a=a';end

% Verify the domain %
if(max(max(u)) > 1 | max(max(v)) > 1 |...
   min(min(u)) < 0 | min(min(v)) < 0)
    error('Values outsite interval [0,1] were found')
end

if(min(min(a)) < 0 )
    error('Negative parameters were found')
end

% Calculate cdf
x = -log(u); y = -log(v);
lxy = log(x) - log(y);
ia = 1./a;
a2 = a./2; 
cdf = exp(-x .* normcdf(ia + a2 .* lxy) - y .* normcdf(ia - a2 .* lxy));
    
%Set exception%
cdf(a==0) = u(a==0) .* v(a==0);
cdf(u==0 & v == 0) = 0;
cdf(u==1 & v == 1) = 1;
    
end % function