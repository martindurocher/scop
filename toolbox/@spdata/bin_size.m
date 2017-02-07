function n = bin_size(obj,varargin)
% SUMMARY : Return the size of the group of lag distances
%
% CALL: obj = obj.bin_size()

if(isempty(varargin))
    n = zeros(1,obj.nbin); 
    for k =1:obj.nbin
        n(k) = size(obj.bins{k},2);
    end
    
else
    n=size(obj.bins{varargin{1}},2);
end