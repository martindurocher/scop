function obj = update(obj,param)
%% Method: update
%
% SUMMARY: Update the parameter of an existing spmodel. 
%
% CALL: out = spmodel.update(param)
%   
% NOTE: see spmodel constructor

switch obj.mdl
    case 'ind'
        error('Independent model does not have parameters')
    case 'normal'
        obj.range = param(1);
        obj.nugget = param(2);
        obj.smooth = param(3);
    case 'chisq'
        obj.a = param(1);
        obj.range = param(2);
        obj.nugget = param(3);
        obj.smooth = param(4);
    case 'student'
        obj.df = param(1);
        obj.range = param(2);
        obj.nugget = param(3);
        obj.smooth = param(4);
    case 't2'
        obj.df = param(1);
        obj.a = param(2);
        obj.range = param(3);
        obj.nugget = param(4);
        obj.smooth = param(5);
    case 'smith'
        obj.cov11 = param(1);
        obj.cov12 = param(2);
        obj.cov22 = param(3);
    otherwise
        error('not implemented')
end
