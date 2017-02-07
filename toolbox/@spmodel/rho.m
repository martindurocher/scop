function out = rho(obj,x)
%% Method: rho
%
% SUMMARY: Return the dependence parameter in respect of the distance. 
%  Except for the Smith model, it correspond to a correlation coefficient
%
% CALL: out = spmodel.rho(x)
%           
% x : Distance where the 
%
switch obj.mdl
    case 'ind'
        out = zeros(size(x));
    case 'smith';
        if(strcmp(obj.opt,'iso'))
            out = 2*sqrt(obj.cov11) ./ x;
        else
           isigma = inv([obj.cov11,obj.cov12;obj.cov12,obj.cov22]);  
           out = 2 ./ sqrt(x*isigma*x');
        end 
        
    case {'normal','chisq','student','t2'}
        out = cor_model(x,obj.link,...
            [obj.range,obj.nugget,obj.smooth]);
        
end 