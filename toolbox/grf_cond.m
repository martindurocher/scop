function [y, sigma,isigma22] = grf_cond(n, coord, new, z, link, ...
                                    param,sigma,isigma22)
                                
m = size(new,1);

% simulate an unconditional grf
if(nargin < 7)
    sigma = cor_model(dist([new',coord']),link,param);
end

uncond = rmvnorm(n,sigma);

%compute simple kriging weights
sigma12 = sigma(1:m,(m+1):end);

if(nargin < 8)
    sigma22 = sigma((m+1):end,(m+1):end);
    isigma22 = pinv(sigma22);
end

w = sigma12*isigma22;

%compute conditional observation
res = repmat(z,1,n)' - uncond(:,(m+1):end);
y = uncond(:,1:m) + res*w';