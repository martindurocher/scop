function pdf = dbhusler(u,a,uselog)
%
%% SUMMARY
%
% Density function of the bivariate Husler-Reiss Copula
%
%% CALL
%
%    pdf = dbhusler_copula(u,a)
%
% u : Uniform sample (Matrix)
% a : Dependence parameter (a>0) 
% uselog : If the log-denstiy should be return
%
%% Author: <martin.durocher@uqtr.ca>

% default parameter %
if(nargin < 3) uselog = 0; end;

v = u(:,2);
u = u(:,1);

[m,n] = size(a);
if(m == 1 & n > 1); a=a';end

% Verify the domain %
if(max(max(u)) >= 1 | max(max(v)) >= 1 |...
   min(min(u)) <= 0 | min(min(v)) <= 0)
    error('Values outsite interval [0,1] were found')
end

x = -log(u); y = -log(v);
lxy = log(x) - log(y);
ia = 1./a;
a2 = a./2; 
mx = ia + a2 .* lxy;
nx = normcdf(mx);
dx = normpdf(mx);
ny = normcdf(ia - a2 .* lxy);

pdf = exp(-x .* nx - y .* ny) ./ (u .* v);
pdf = pdf .* ( (ny .* nx) + (a2 ./ y) .* dx);

% Set exception %
pdf(a==0) = 1;

if(uselog)
    pdf = log(pdf);
end

end % function