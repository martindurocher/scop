function [theta, nllik] = fit_pairwise(h,u,sp)
%%
%% SUMMARY:
%
%    Calculate by pairwise likelihood the range parameter for 
%          a copula model
%
%% CALL
%
%    [theta, nllik] = range_fit_pairwise(h,u,sp)
%  
%  theta : Range of the link function (Numeric)
%  nllik : Negative log likelihood
%  h     : Distance (Matrix n x 1)
%  u     : ranks of the pairs (Matrix n x 2)
%  sp    : a spmodel object
%
%% NOTE
%
% see cor_model.m for description of the correlation model
% see sim_copula_multi for a description of the copula
%
% The variable mdl specify the model. By default, only the 
% range parameter is estimated. For instance, mdl = 'normal'
% estimate the range parameter of a normal copula. The suffix
% '_n','_s', '_b' can be added to estimate also the nugget, 
% the smooth or both parameter. For instance mdl = 'normal_n' will
% estimate the normal copula with a range and a nugget effect.
% idem for chisq. For Smith model there is the isotropic field 
% ('smith_i') or the default anysotropic field ('smith').
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>

mdl = [sp.mdl,'_',sp.opt];
link = sp.link;
param = sp.theta();

switch link
    case 'pow'
        p3 = 1;
    otherwise;
        p3 = 0;
end

switch mdl

    case 'normal_r'
        param(1) = log(param(1));
        fobj = @(p) npllik_bnorm([exp(p),param(2),param(3)],...
                link,u,h);
        [theta, nllik] = fminsearch(fobj,param(1),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        
        theta = [exp(theta),param(2),param(3)];
        
    case 'normal_rn'
        param(1) = log(param(1));
        param(2) = logit(param(2));
                
        fobj = @(p) npllik_bnorm([exp(p(1)),expit(p(2)),param(3)],...
                link,u,h);
            
        [theta, nllik] = fminsearch(fobj,[param(1),0],...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        theta = [exp(theta(1)),expit(theta(2)),param(3)];
        
    case 'normal_rs'
        if(not(p3)); 
            error('Must be a three parameter link function')
        end
        param(1) = log(param(1));        
        
        fobj = @(p) npllik_bnorm([exp(p(1)),param(2),p(2)],...
                link,u,h);

        [theta, nllik] = fminsearch(fobj,param([1,3]),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        
        theta = [exp(theta(1)),param(2),theta(2)];
                  
    case 'chisq_r'
        param(2) = log(param(2));
        
        fobj = @(p) npllik_bchisq([exp(p),param(3),param(4)],...
                link,param(1),u,h,sp.tab_ha);  
            
        [theta, nllik] = fminsearch(fobj,param(2),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        theta = [param(1),exp(theta),param(3),param(4)];
        
    case 'chisq_rn'
        param(2) = log(param(2));
        param(3) = logit(param(3));
                
        fobj = @(p) npllik_bchisq([exp(p(1)),...
                expit(p(2)) , param(4)],...
                link,param(1),u,h,sp.tab_ha);
            
        [theta, nllik] = fminsearch(fobj,param(2:3),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        theta = [param(1),exp(theta(1)),expit(theta(2)),param(4)];
        
    case 'chisq_rs'
        if(not(p3)); 
            error('Must be a three parameter link function')
        end
        param(2) = log(param(2));
               
        fobj = @(p) npllik_bchisq([exp(p(1)),param(3),p(2)],...
                link,a,u,h,sp.tab_ha);
            
        [theta, nllik] = fminsearch(fobj,param([1,3]),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        
        theta = [param(1),exp(theta(1)),param(3),theta(2)];
       
    case 'chisq_ar'    
        param(1) = log(param(1));
        param(2) = log(param(2));
                       
        fobj = @(p) npllik_bchisq([exp(p(2)),param(3),param(4)],...
                link,exp(p(1)),u,h,[]);
            
        [theta, nllik] = fminsearch(fobj,param(1:2),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        
        theta = [exp(theta(1)),exp(theta(2)),param(3:4)];

    case 'student_r'
        param(2) = log(param(2)); 
        
        u = tinv(u,param(1));
        
        fobj = @(p) npllik_bt([exp(p),param(3),param(4)],...
            link,param(1),u,h);
        
        [theta, nllik] = fminsearch(fobj,param(2),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        
        theta = [param(1),exp(theta(1)),param(3:4)];
        
    case 'student_rn'
        param(2) = log(param(2));
        param(3) = logit(max(.01,param(3)));
        
        u = tinv(u,param(1));
        
        fobj = @(p) npllik_bt([exp(p(1)),expit(p(2)),param(4)],...
            link,param(1),u,h);
        
        [theta, nllik] = fminsearch(fobj,param(2:3),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        
        theta = [param(1),exp(theta(1)),expit(theta(2)),param(4)];
        
    case 't2_r'
        %[nu,a,range,nugget,smooth]
        param(3) = log(param(3)); 
        
        up = (1 + u(:,1))/2; 
        un = (1 - u(:,1))/2; 
        v  = (1 + u(:,2))/2;
        u = tinv([up,un,v],param(1));
        
        fobj = @(p) npllik_bt2([exp(p),param(4),param(5)],...
            link,param(1),param(2),u,h);
        
        [theta, nllik] = fminsearch(fobj,param(3),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        
        theta = [param(1:2),exp(theta(1)),param(4:5)];

    case 't2_rn'
        error('Not yet')
        
    case 't2_dr'
        error('Not yet')
        
    case 't2_ar'
        error('Not yet')
        
    case 'smith_iso'
        param(1) = 1./log(param(1));
        fobj = @(p) npllik_bsmith_i(exp(1./p),u,h);
        [theta, nllik] = fminsearch(fobj,param(1),...
            optimset('MaxFunEvals', 5000, 'MaxIter', 5000));
        
        theta = [exp(1./theta),0,exp(1./theta)];
        
    case 'smith_ani'
        error('not implemented')
            
end

end

function L = npllik_bnorm(theta_link,link,u,h)
    rho = cor_model(h,link,theta_link);
    L = -sum(copula.dbnorm(u,rho,1));
end    

function L = npllik_bchisq(theta_link,link,a,u,h,tab)
    rho = cor_model(h,link,theta_link);
    L = -sum(copula.dbchisq(u,rho,a,1,tab));    
end
       
function L = npllik_bsmith_i(theta_link,u,h)
    rho = 2*sqrt(theta_link)./h;
    L = -sum(copula.dbhusler(u,rho,1));
end

function L = npllik_bt(theta_link,link,nu,u,h)
    rho = cor_model(h,link,theta_link);
    L = -sum(copula.dbt(u,rho,nu,1,0));
end 

function L = npllik_bt2(theta_link,link,nu,a,u,h)
    rho = cor_model(h,link,theta_link);
    L = -sum(copula.dbt2(u,rho,nu,a,1,0));
end 