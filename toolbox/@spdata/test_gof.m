function [out,sp0] = test_gof(obj, type, nboot, method, h0)

% default parameters
if(nargin < 2); type  = 'xchk'; end
if(nargin < 3); nboot = 0; end
if(nargin < 4); method = 'pairwise'; end
if(nargin < 5); h0 = 'composite'; end
    
% default parameters 
if(isempty(obj.nbin))
   error('Breaks must be specified') 
end

%Organizing arguments for various methods
switch h0
    
    %in case of independence
    case {'composite'}
        theo_type = 1;
        
    case {'table'}
        h0 = 'composite';
        theo_type = 3;
        
    otherwise
        error('Wrong hypothesis H0')
end

% Declare variables%
sp0 = [];
stat = struct('cvm',[],'cvmk',[],'ks',[],'ksk', [], ...
    'xch', [], 'xchk', [], 'cvmkt', [], 'xchkt', []);
boot = struct('cvm',[],'cvmk',[],'ks',[],'ksk', [], ...
    'xch', [], 'xchk', [], 'cvmkt', [], 'xchkt', []);
do_stat = struct('cvm',0,'cvmk',0,'ks',0,'ksk',0,...
    'xch',0,'xchk',0, 'cvmkt', 0, 'xchkt', 0);

%Determine which gof stats to calculate %  
switch type
    case 'cvm'
        do_stat.cvm = 1;
        stat.cvm = zeros(1,obj.nbin);
        boot.cvm = zeros(nboot,obj.nbin);
        
    case 'cvmk'
        do_stat.cvmk = 1;
        stat.cvmk = zeros(1,obj.nbin);
        boot.cvmk = zeros(nboot,obj.nbin);
        
    case 'ks'
        do_stat.ks = 1;
        stat.ks = zeros(1,obj.nbin);
        boot.ks = zeros(nboot,obj.nbin);
    case 'ksk'
        do_stat.ksk = 1;
        stat.ksk = zeros(1,obj.nbin);
        boot.ksk = zeros(nboot,obj.nbin);
        
    case 'xch'
        do_stat.xch = 1;
        stat.xch = zeros(1,obj.nbin);
        boot.xch = zeros(nboot,obj.nbin);
        
    case 'xchk'
        do_stat.xchk = 1;
        stat.xchk = zeros(1,obj.nbin);
        boot.xchk = zeros(nboot,obj.nbin);                   
        
    case 'cvm2'
        do_stat.cvm = 1;
        stat.cvm = zeros(1,obj.nbin);
        boot.cvm = zeros(nboot,obj.nbin);
        
        do_stat.cvmk = 1;
        stat.cvmk = zeros(1,obj.nbin);
        boot.cvmk = zeros(nboot,obj.nbin); 
        
    case 'ks2'
        do_stat.ks = 1;
        stat.ks = zeros(1,obj.nbin);
        boot.ks = zeros(nboot,obj.nbin);
        
        do_stat.ksk = 1;
        stat.ksk = zeros(1,obj.nbin);
        boot.ksk = zeros(nboot,obj.nbin);
        
    case 'xch2'
        do_stat.xch = 1;
        stat.xch = zeros(1,obj.nbin);
        boot.xch = zeros(nboot,obj.nbin);
        
        do_stat.xchk = 1;
        stat.xchk = zeros(1,obj.nbin);
        boot.xchk = zeros(nboot,obj.nbin);
        
    case 'cvmx'
        do_stat.xch = 1;
        stat.xch = zeros(1,obj.nbin);
        boot.xch = zeros(nboot,obj.nbin);
        
        do_stat.xchk = 1;
        stat.xchk = zeros(1,obj.nbin);
        boot.xchk = zeros(nboot,obj.nbin);
        
        do_stat.cvm = 1;
        stat.cvm = zeros(1,obj.nbin);
        boot.cvm = zeros(nboot,obj.nbin);
        
        do_stat.cvmk = 1;
        stat.cvmk = zeros(1,obj.nbin);
        boot.cvmk = zeros(nboot,obj.nbin);
        
    case 'all'
        do_stat.cvm = 1;
        stat.cvm = zeros(1,obj.nbin);
        boot.cvm = zeros(nboot,obj.nbin);
        
        do_stat.cvmk = 1;
        stat.cvmk = zeros(1,obj.nbin);
        boot.cvmk = zeros(nboot,obj.nbin);
        
        do_stat.ks = 1;
        stat.ks = zeros(1,obj.nbin);
        boot.ks = zeros(nboot,obj.nbin);
        
        do_stat.ksk = 1;
        stat.ksk = zeros(1,obj.nbin);
        boot.ksk = zeros(nboot,obj.nbin);
        
        do_stat.xch = 1;
        stat.xch = zeros(1,obj.nbin);
        boot.xch = zeros(nboot,obj.nbin);

        do_stat.xchk = 1;
        stat.xchk = zeros(1,obj.nbin);
        boot.xchk = zeros(nboot,obj.nbin);
        
    otherwise
        disp('No GOF selected. Used as bootstrap procedure');
end

% Deterimine the hypothesis h0 and estimate parameter %
switch h0
    case 'table'
        sp0 = obj.sp;  
    case 'composite'
        sp0 = obj.estimate(method);
end

% compute the parameter of each bivariate copula
rho = sp0.rho(obj.mids);

% Calculate test statistics
if(~strcmp(h0,'table')) % not necessary for building a table
for k = 1:obj.nbin 
    
    %sample from H0 process to approximate theoritical copula
    theo = sp0.sim_biv(obj.GOF_N_APPROX,rho(k));
    
    % compute the tests statistics
    if(do_stat.cvm)
        stat.cvm(k) = copula.distance(obj.get_pairs(k),...
            theo,'cvm','emp',obj.HISTSIZE);
    end
    if(do_stat.cvmk)
        stat.cvmk(k) = copula.distance(obj.get_pairs(k),...
            theo,'cvm','kendall');
    end
    if(do_stat.ks)
        stat.ks(k) = copula.distance(obj.get_pairs(k),...
            theo,'ks','hist',obj.HISTSIZE);
    end    
    if(do_stat.ksk)
        stat.ksk(k) = copula.distance(obj.get_pairs(k),...
            theo,'ks','kendall');
    end 
    if(do_stat.xch)
        stat.xch(k) = copula.distance(obj.get_pairs(k),...
            theo,'xch','emp',obj.HISTSIZE);
    end    
    if(do_stat.xchk)
        stat.xchk(k) = copula.distance(obj.get_pairs(k),...
            theo,'xch','kendall');
    end
    
end %for
end %if not a table

% If a table is used later for calculating the critical value
if(nboot == 0)
   out = gof_copula(stat,boot,type); 
   return 
end
 
%otherwise parametric bootstrap start here
for ii=1:nboot
    
    % for monitoring: text bar or random point%
    if(obj.PROGRESS)
        monitor(ii,nboot)             
    end
            
    % bootstrap sample %
    switch sp0.mdl
        case 'smith'
            warning('Not tested')
            bsample = sp0.sim_coord(1,obj.coord)';
        otherwise
            bsample = sp0.sim_dist(1,obj.distance)';     
    end
    
    % calculate the ranks if required
    if(strcmp(obj.type,'ranks'))
        bsample = rankit(bsample);
    end
    
    boot_data = pairwise(bsample);
    
    % estimate parameters for the bootrap sample
    switch method
        case 'full'
            spb = obj.estimate('full',bsample);
        case 'pairwise'
            spb = obj.estimate('pairwise',boot_data); 
    end
    
    rho = spb.rho(obj.mids); 
    
    % Calculate test the statistics %
    for k = 1:obj.nbin        
        
        theo = spb.sim_biv(obj.GOF_N_APPROX,rho(k));       
    
        %compute the test statistics
        if(do_stat.cvm)
            boot.cvm(ii,k) = copula.distance(...
                boot_data(obj.bins{k},:),...
                theo,'cvm','emp',obj.HISTSIZE);
        end
        if(do_stat.cvmk)
            boot.cvmk(ii,k) = copula.distance(...
                boot_data(obj.bins{k},:),...
                theo,'cvm','kendall');
        end
        if(do_stat.ks)
            boot.ks(ii,k) = copula.distance(...
                boot_data(obj.bins{k},:),...
                theo,'ks','hist',obj.HISTSIZE);
        end    
        if(do_stat.ksk)
            boot.ksk(ii,k) = copula.distance(...
                boot_data(obj.bins{k},:),...
            theo,'ks','kendall');
        end
        if(do_stat.xch)
            boot.xch(ii,k) = copula.distance(...
                boot_data(obj.bins{k},:),...
                theo,'xch','emp',obj.HISTSIZE);
        end    
        if(do_stat.xchk)
            boot.xchk(ii,k) = copula.distance(...
                boot_data(obj.bins{k},:),...
            theo,'xch','kendall');
        end 
        
    end
                    
end

% Output bootstrap objects
out = gof_copula(stat,boot,type);

 