function obj = set_breaks(obj, breaks, type)
%
% SUMMARY : set the breaks that define the group of lag distances
%
% CALL: obj = obj.set_breaks(breaks,type)
%
%       breaks : inner delimiter (Numeric 1 x m)
%     
%       type : 'mid' = middel of the intervals)
%              'med' = median of the paired distances
%              'avg' = average of the paired distances
%

obj.breaks = breaks;
obj.bins = binning(obj.pdistance,breaks);
obj.nbin = size(obj.bins,1);

% Verify that all bin are not empty
for(ii = 1:obj.nbin)
   if(size(obj.bins{ii}, 2) == 0)
       error('An empty bin was created');
   end
end

rep = zeros(1,obj.nbin);
switch type
    case 'mid'
        rep = (breaks(2:(obj.nbin+1)) + obj.breaks(1:obj.nbin))/2;
    case 'med'
        for k = 1:obj.nbin
            rep(k) = median(obj.pdistance(obj.bins{k}));
        end             
    case 'avg'
        for k = 1:obj.nbin
            rep(k) = mean(obj.pdistance(obj.bins{k}));
        end
end
obj.mids = rep;
