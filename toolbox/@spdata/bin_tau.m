function tau = bin_tau(obj)
% SUMMARY : Return kendall tau of the group of lag distances
%
% CALL: obj = obj.bin_tau()

tau = zeros(1,obj.nbin);
for k=1:obj.nbin
    tau(k) = copula.tau(obj.get_pairs(k));
end       
