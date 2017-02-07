function y = eff_range(obj,pars,span,D)
%
%% SUMMARY
%
% Function to calculate the effective range of a spatial copula model
%
%% CALL
%
%          y = eff_range(pars,span, D)
%
% span: interval containing the effective range
% link: Link function
% pars: Parameter of the link function (range,nugget, smooth)
% mdl : Copula model
% a:    Copula specific parameter 
% D: value of tau at the effective range (default = 0.2)

%% SEE ALSO
%
% cor_model.m, spmodel.theta_tau.m
%
%% AUTHOR: martin du rocher <martin.durocher@uqtr.ca>
%

%Default value
if(nargin < 3); span = [.01,obj.range]; end
if(nargin < 4); D = .2; end

if(nargin <2)
    pars = obj.theta('c');
end

switch obj.mdl
    case {'normal','student'}  
        fobj = @(x) copula.theta_tau(...
                cor_model(x,obj.link,pars),obj.mdl)-D;
            
    case 'chisq'    
        fobj = @(x) copula.theta_tau(...
            cor_model(x,obj.link,pars),obj.mdl,obj.a)-D;
end

y = fzero(fobj,span);