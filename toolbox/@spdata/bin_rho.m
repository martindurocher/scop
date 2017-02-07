function srho = bin_rho(obj)
% SUMMARY : Return Spearman rho of the group of lag distances
%
% CALL: obj = obj.bin_rho()

     srho = zeros(1,obj.nbin);
     for k=1:obj.nbin
         tmp = corr(obj.pairs(obj.bins{k},:));
         srho(k) = tmp(2);
     end       
end