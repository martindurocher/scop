function y = get_pairs(obj,k)

if(isempty(obj.nbin))
   error('Breaks must be specified') 
end

if(k==0)
    y = obj.pairs;
else
    y = obj.pairs(obj.bins{k},:);
end
