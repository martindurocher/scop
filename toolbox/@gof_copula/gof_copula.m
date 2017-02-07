classdef gof_copula
%%
%% OBJECTIVE 
%
%  Produce the summary of the output coming form spdata.test_gof()
%     and spdata.test_ind()
%
%% CALL
%
% [stat,pval,test] = summary(obj,type)
%
% stat : Test statistics for different groups of lags and 
%        global statistics (Vector)
% pval : P-values of the satistics (Vector)
% test : Conclusion of the test (Vector)
% type : Which statistics to summarize (String)
%
%% NOTE
% 
% Avalaible type:

% Test based on empirical copula: 
%       Cramer Von Mises = 'cvm',
%       'cvm' + symmetry = 'xch'
%       Kolmogorov-Smirnov = 'ks'

% Test based on Kendall transformation
%       Cramer Von Mises = 'cvmk',
%       CVM + symmetry = 'xchk'
%       Kolmogorov-Smirnov = 'ksk'
%
% Combination of type 
%  'cvm2' = 'cvm'+ 'cvmk' (idem for 'xch' and 'ks')
%  'all' = all tests are performed
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>

properties
    cvm  % Sample CVM
    cvmk  % sample CVM with kendall transform
    cvmkt 
    ks  % sample KS
    ksk  % sample KS with kendall transform
    xch
    xchk
    xchkt
    
    % Bootstrap list of the respective statistics
    Bcvm 
    Bcvmk
    Bcvmkt
    Bks
    Bksk
    Bxch
    Bxchk
    Bxchkt
    
    nbin
    nboot
    
    % Type of test statistics
    type
    
    % For external bootstrap table
    table
    id_size
end
  
methods
  
function obj = gof_copula(stat,boot,type)
    obj.cvm = stat.cvm;
    obj.cvmk = stat.cvmk;
    obj.cvmkt = stat.cvmkt;
    
    obj.ks = stat.ks;
    obj.ksk = stat.ksk;
    
    obj.xch = stat.xch;
    obj.xchk = stat.xchk;
    obj.xchkt = stat.xchkt;
    
    obj.Bcvm = boot.cvm;
    obj.Bcvmk = boot.cvmk;
    obj.Bcvmkt = boot.cvmkt;
    
    obj.Bks = boot.ks;
    obj.Bksk = boot.ksk;
    
    obj.Bxch = boot.xch;
    obj.Bxchk = boot.xchk;
    obj.Bxchkt = boot.xchkt;
    
    obj.nbin = max([size(obj.cvm,2),size(obj.cvmk,2),...
        size(obj.ks,2),size(obj.ksk,2),...
        size(obj.xch,2),size(obj.xchk,2),...
        size(obj.cvmkt,2),size(obj.xchkt,2)]);
    
    obj.nboot = max([size(obj.Bcvm,1),size(obj.Bcvmk,1),...
        size(obj.Bks,1),size(obj.Bksk,1),...
        size(obj.Bxch,1),size(obj.Bxchk,1),...
        size(obj.Bcvmkt,1),size(obj.Bxchkt,1)]);
   
    obj.type = type;
    
end
    
end %methods

end %class    