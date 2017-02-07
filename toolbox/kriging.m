function [zpred,vpred,w] = kriging(z,coord,new,link,param,opt)

if(nargin < 6); opt = 'simple'; end;

m = size(new,1);
n = size(coord,1);

% transpose z if necessary
if(size(z) ~= [n,1])
   error('wrong dimension for Z')
end

sigma = cor_model(dist([new',coord']),link,param);

sigma11 = sigma(1:m,1:m);
sigma12 = sigma(1:m,(m+1):end);
sigma22 = sigma((m+1):end,(m+1):end);
isigma22 = inv(sigma22);
w = sigma12*isigma22;

switch opt
    case {'s','simple'}
        zpred = (w*z)';
        vpred = sigma11 - w*sigma12';
    case {'o','ordinaire'}
        %one = ones(n,1);
        %sig1 = one'*isigma22;
        %isig = (sig1*one);
        %sig2 = sig1*sigma12';
        
        %mu =  sig1*z / isig;
        %zpred = mu + w * (z - one*mu);
        
        %lambda = 1 - sig2 / isig;
        
        %spred = sigma11 - w*sigma12';
end

