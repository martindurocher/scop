function u = sim_cond(obj,n,new,m)

if(nargin < 4) m = []; end;

nnew = size(new,1);
param = obj.sp.theta();

switch obj.sp.mdl
    case 'normal'
        z0 = norminv(obj.ranks);
        zcond = grf_cond(n, obj.coord, new, z0,obj.sp.link, param);
        u = normcdf(zcond);
        
    case 'chisq'   
        if(param(1) == 0)
            z0 = norminv( (obj.ranks + 1) / 2);
            zcond = grf_cond(n,obj.coord,new,z0,...
                obj.sp.link,param(2:4));
            u = 2 * normcdf( abs(zcond+param(1)) ) - 1;
        else    
            if(~isempty(obj.sp.tab_ha))
                ha = interp1(obj.sp.tab_ha(:,1),...
                    obj.sp.tab_ha(:,2),...
                    obj.ranks,'spline');
            else
                ha = copula.chisq_ha(obj.ranks,param(1));
            end 
            
            z0 = -norminv(1-ha);
                      
            zcond=  grf_cond(n,obj.coord,new,...
               z0, obj.sp.link,param(2:4));
            
            y = abs(zcond+param(1));
            u = normcdf(y-param(1)) + normcdf(y+param(1)) - 1;
            
        end
        
    case 'student'
        w = sqrt(param(1)./chi2rnd(param(1),1,n));
        u = zeros(n,nnew);
                    
        for(kk = 1 :n)
           z0 = tinv(obj.ranks,param(1))/w(kk);
           
           if(kk == 1)
                [zcond,sig,isig22] =  grf_cond(1,obj.coord,new,...
                                    z0, obj.sp.link,param(2:4));
           else
               zcond =  grf_cond(1,obj.coord,new, z0, ...
                   obj.sp.link,param(2:4),sig,isig22);
           end
           
           u(kk,:) = tcdf( w(kk)*zcond, param(1));
        end
    
    otherwise
        error('Conditional simulation is not avalaible for that copula')
end