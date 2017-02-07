function out = theta(obj,opt)
%% Method: theta
%
% SUMMARY: Return the a vector of the parameter of the model. 
%
% CALL: out = spmodel.theta()
%           
% EXAMPLE: 
%  sp0 = spmodel('normal','exp',[.3,0,1],'r')
%  [range0,nugget0,smooth0] = sp0.theta()
%
% NOTE: see spmodel constructor
% opt = d : return all parameter
% opt = a : return copula-specific parameter
% opt = c : return correlation parameter
% opt = n : return the number of parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin <2); opt = 'd'; end

model = [obj.mdl,'_',opt];

switch model
    case {'ind_d','ind_c','ind_a'}
        error('Independent model does not have parameters');
    case {'normal_d','normal_c'}
        out = [obj.range,obj.nugget,obj.smooth];
    case 'normal_a'
        out = [];
    case 'chisq_d'
        out = [obj.a,obj.range,obj.nugget,obj.smooth];
    case 'chisq_c'
        out = [obj.range,obj.nugget,obj.smooth];
    case 'chisq_s'
        out = [obj.a];
    case 'smith_d'
         out = [obj.cov11,obj.cov12,obj.cov22];
    case 'student_d'
         out = [obj.df,obj.range,obj.nugget,obj.smooth];
    case 'student_c'
         out = [obj.range,obj.nugget,obj.smooth];
    case 'student_a'
         out = [obj.df];
    case 't2_d'
         out = [obj.df,obj.a,obj.range,obj.nugget,obj.smooth];
    case 'normal_n'
        out = 3;
    case {'chisq_n','student_n'}
        out = 4;
    case 't2'
        out = 5;
    otherwise
        error('not implemented')
end
