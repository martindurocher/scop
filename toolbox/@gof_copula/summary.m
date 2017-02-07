function [stat,pval,test] = summary(obj, type, alpha, tab)

%default value
if(nargin < 3); alpha = .05; end
if(nargin < 2); type = obj.type; end
do_tab = (nargin == 4);

% select desired statistics
switch type
    case 'cvm'
        stat = obj.cvm;
        boot = obj.Bcvm;
    case 'cvmk'
        stat = obj.cvmk;
        boot = obj.Bcvmk;
    
    case 'cvmkt'
        stat = obj.cvmkt;
        boot = obj.Bcvmkt;
    
    case 'ks'
        stat = obj.ks;
        boot = obj.Bks;
    case 'ksk'
        stat = obj.ksk;
        boot = obj.Bksk;
        
    case 'xch'
        stat = obj.xch;
        boot = obj.Bxch;
    case 'xchk'
        stat = obj.xchk;
        boot = obj.Bxchk;
    case 'xchkt'
        stat = obj.xchkt;
        boot = obj.Bxchkt; 
end


if(do_tab)
    stat = [stat, cumsum(stat)];
    test = (stat < tab);
    pval = [];
else
    boot = [boot,cumsum(boot')'];      
    stat = [stat, cumsum(stat)];
        
    nk = max(size(stat));
    pval = zeros(1,nk);
    
    for k=1:nk       
        pval(k) = mean(boot(:,k)>stat(k));
    end 

    test = (pval > alpha); 
end




