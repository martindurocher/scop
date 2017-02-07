function [obj, bootp]= fit(obj,varargin)

%%%%%%%%% HEADER %%%%%%%%%%%%%%%%%%%%%%%
method = 'pairwise';
nboot = 0;
htol = obj.HTOL;
bootp = 0;

for(ii = 1:numel(varargin))
switch varargin{ii}
    case {'pairwise','full'}
        method = varargin{ii};
    case {'nboot','bootstrap'}
        nboot = varargin{ii+1};
    case 'htol'
        htol = varargin{ii+1};
end
end

obj.HTOL = htol;
fit = obj.estimate(method);
obj.sp = fit;

% bootstrap parameter if necessary
if(nboot>0)
    bootp = zeros(nboot,fit.theta('n'));
    sample = fit.sim_coord(nboot,obj.coord);

    for(kk = 1:nboot)
        monitor(kk,nboot);
        x = sample(kk,:)';
        switch method
            case 'full'
                boot_fit = obj.estimate('full',x);
            case 'pairwise'
                boot_fit = obj.estimate('pairwise',pairwise(x)); 
        end
        bootp(kk,:) = boot_fit.theta();
   
    end
end