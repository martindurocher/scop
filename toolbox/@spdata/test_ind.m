function out = test_ind(obj, type, nboot)

% default parameters
if(nargin < 2); type  = 'cvmk'; end
if(nargin < 3); nboot = 0; end

if(isempty(obj.nbin))
   error('Breaks must be specified') 
end

% Declare variables%
stat = struct('cvm',[],'cvmk',[],'cvmkt',[],'ks',[],'ksk', [], ...
    'xch', [], 'xchk', [], 'xchkt',[]);
boot = struct('cvm',[],'cvmk',[],'cvmkt',[],'ks',[],'ksk', [], ...
    'xch', [], 'xchk', [], 'xchkt', []);
do_stat = struct('cvm',0,'cvmk',0,'cvmkt',0,'ks',0,'ksk',0,...
    'xch',0,'xchk',0,'xchkt',0);

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
        
    case 'cvmkt'
        do_stat.cvmkt = 1;
        stat.cvmkt = zeros(1,obj.nbin);
        boot.cvmkt = zeros(nboot,obj.nbin);
        
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
    
    case 'xchkt'
        do_stat.xchkt = 1;
        stat.xchkt = zeros(1,obj.nbin);
        boot.xchkt = zeros(nboot,obj.nbin);  
        
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

for k = 1:obj.nbin   
    if(do_stat.cvm)
        stat.cvm(k) = copula.distance_ind(...
            obj.get_pairs(k),'cvm');
    end
    if(do_stat.cvmk)
        stat.cvmk(k) = copula.distance_ind(...
            obj.get_pairs(k),'cvmk');
    end
    if(do_stat.cvmkt)
        stat.cvmkt(k) = copula.distance_ind(...
            obj.get_pairs(k),'cvmkt');
    end
    if(do_stat.ks)
        stat.ks(k) = copula.distance_ind(...
            obj.get_pairs(k),'ks');
    end    
    if(do_stat.ksk)
        stat.ksk(k) = copula.distance_ind(...
            obj.get_pairs(k),'ksk');
    end 
    if(do_stat.xch)
        stat.xch(k) = copula.distance_ind(...
            obj.get_pairs(k),'xch');
    end    
    if(do_stat.xchk)
        stat.xchk(k) = copula.distance_ind(...
            obj.get_pairs(k),'xchk');
    end
    if(do_stat.xchkt)
        stat.xchkt(k) = copula.distance_ind(...
            obj.get_pairs(k),'xchkt');
    end
    
end %for

% If a table is used for calculating the critical value
if(nboot == 0)
   out = gof_copula(stat,boot,type); 
   return 
end

for ii=1:nboot
    
    switch obj.PROGRESS
        case 2
            pstep = max(1,round(nboot/20));
            if(mod(ii,5*pstep)==0); fprintf('|');
            elseif(mod(ii,pstep)==0); fprintf('.'); end
        case 1
            pause(0.00001)
            if(rand(1)<(.005/nboot)); fprintf(1,'\b.\n'); end;                
    end
    
    boot_data = pairwise(rand(obj.nsite,1));
    
    for k = 1:obj.nbin
        bdata = boot_data(obj.bins{k},:);
            
        if(do_stat.cvm)
            boot.cvm(ii,k) = copula.distance_ind(bdata,'cvm');
        end
        if(do_stat.cvmk)
            boot.cvmk(ii,k) = copula.distance_ind(bdata,'cvmk');
        end
        if(do_stat.cvmkt)
            boot.cvmk(ii,k) = copula.distance_ind(bdata,'cvmk');
        end
        if(do_stat.ks)
            boot.ks(ii,k) = copula.distance_ind(bdata,'ks');
        end    
        if(do_stat.ksk)
            boot.ksk(ii,k) = copula.distance_ind(bdata,'ksk');
        end
        if(do_stat.xch)
            boot.xch(ii,k) = copula.distance_ind(bdata,'xch');
        end    
        if(do_stat.xchk)
            boot.xchk(ii,k) = copula.distance_ind(bdata,'xchk');
        end
        if(do_stat.xchkt)
            boot.xchkt(ii,k) = copula.distance_ind(bdata,'xchkt');
        end 
        
    end
    
end

if(obj.PROGRESS == 2); fprintf('\n'); end;

% Output bootstrap objects
out = gof_copula(stat,boot,type);


