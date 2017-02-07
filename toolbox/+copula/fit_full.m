function [theta, nllik] = fit_full(h,u,sp)
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
        
        x = norminv(u');
        
        fobj = @(p) nllik_mvnorm([exp(p),param(2),param(3)],...
                link,x,h);
            
        [theta, nllik] = fminsearch(fobj,param(1));
        theta = [exp(theta),param(2),param(3)];
        
    case 'normal_rn'
        param(1) = log(param(1));
        param(2) = logit(param(2));
        
        x = norminv(u');
        
        fobj = @(p) nllik_mvnorm([exp(p(1)),expit(p(2)),param(3)],...
                link,x,h);
            
        [theta, nllik] = fminsearch(fobj,[param(1),0]);
        theta = [exp(theta(1)),expit(theta(2)),param(3)];
        
    case 'normal_rs'
        if(not(p3)); 
            error('Must be a three parameter link function')
        end
        param(1) = log(param(1)); 
        
        x = norminv(u');
        
        fobj = @(p) nllik_mvnorm([exp(p(1)),param(2),p(2)],...
                link,x,h);

        [theta, nllik] = fminsearch(fobj,param([1,3]));
        
        theta = [exp(theta(1)),param(2),theta(2)];
                  
    case 'student_r'
        param(2) = log(param(2)); 
        
        x = tinv(u',param(1));
        
        fobj = @(p) nllik_mvt([exp(p),param(3),param(4)],...
            link,param(1),x,h);
        
        [theta, nllik] = fminsearch(fobj,param(2));
        theta = [param(1),exp(theta(1)),param(3:4)];
        
    case 'student_rn'
        param(2) = log(param(2));
        param(3) = logit(max(.01,param(3)));
        
        x = tinv(u',param(1));
        
        fobj = @(p) nllik_mvt([exp(p(1)),expit(p(2)),param(4)],...
            link,param(1),x,h);
        
        [theta, nllik] = fminsearch(fobj,param(2:3));
        theta = [param(1),exp(theta(1)),expit(theta(2)),param(4)];
            
end

end

function L = nllik_mvnorm(theta_link,link,x,h)
    sigma = cor_model(h,link,theta_link);
    L = -log(mvnpdf(x,zeros(1,size(sigma,1)),sigma));
end    

function L = nllik_mvt(theta_link,link,nu,x,h)
    sigma = cor_model(h,link,theta_link);
    L = -log(mvtpdf(x,sigma,nu));
end 
