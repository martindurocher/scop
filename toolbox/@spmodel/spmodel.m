classdef spmodel
%
%% Summary
%
%  Object containing the functionality of a spatial model using
% the copula-based semiparametric approach.
%
%% List of properties
% * range : range parameter
% * nugget: nugget effect
% * smooth: smooth parameter
% * cov11 : variance 1th variate (smith model) 
% * cov12 : covariance coefficient (smith model) 
% * cov22 : variance 2nd variate (smith model)
% * mdl   : Model name
% * a     : non-centrality parameter
% * df    : degree of freedom (scalar)
% * opt   : Option (char)
% * link  : link function (char)
% * nllik : negative log pairwise likelihood 
% 
%% Method: spmodel
% 
% SUMMARY: constructor
% 
% CALL: obj = spmodel(mdl, link, param, opt)
% 
% * mdl  : Model name
% * link : Name of the link function
% * param: Vector of initial parameters (in respect of mdl)
% * opt  : Additional option (char)
%
% NOTE: 
%
% The model names and parameter are : 
% 'normal'     , [range,nugget,smooth]
% 'chi-squared', [a,range,nugget,smooth]
% 'student'    , [df,range,nugget,smooth]
% 't2'         , [a,df,range,nugget,smooth]
% 'smith'      , [cov11,cov12,cov22]
%
% For mor details see function: cor_model.m
%
% The optional methode indicate which parameter are fixed and 
% variable. For more details see the method estimate
% of the spmodel object.
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>

properties
    range;
    nugget;
    smooth;
    cov11;
    cov12;
    cov22;
    mdl;
    a;
    df;
    opt;
    link;
    nllik=0;
    nllik_type='';
    tab_ha=[];
    
end %properties

methods
    
function obj = spmodel(mdl, link, param, opt)

    %default value
    if(nargin == 1);  param = []; end;
    if(nargin < 4); opt = 'r'; end;
    
    psize = size(param,2);
    
    switch mdl
        case 'ind'
            obj.mdl = mdl;
            
        case 'smith';
                          
            obj.mdl = mdl;
            obj.link = link;
            obj.opt = link;
            
            % validate parameter
            if(param(1) < 1e-6)
                error('Invalid parameter')
            else
                obj.cov11 = param(1);
            end
                        
            switch link
                case 'iso'
                    if(psize~=1)
                        error('Wrong number of parameter')
                    end
                    
                    obj.cov12 = 0;
                    obj.cov22 = param(1);
                case 'ani'
                    if(size(param,2)~=1)
                        error('Wrong number of parameter')
                    end
                    
                    % validate parameter
                    if(param(3) < 1e-6)
                        error('Invalid parameter')
                    elseif(abs(param(2))> 1 )
                        error('Invalid parameter')
                    else    
                        obj.cov12 = param(2);
                        obj.cov22 = param(3);
                    end
            end
                        
        case {'normal'}
                  
            if(psize == 1)
                param = [param,0,1];
            elseif(psize == 2)
                param = [param,1];
            elseif(size(param,2)>3)
                error('Wrong number of parameter')
            end
                    
            obj.mdl = mdl;
            obj.a = 0;
            obj.range = param(1);
            obj.nugget = param(2);
            obj.smooth = param(3);
            obj.opt = opt;
            obj.link = link;
            
        case {'chisq'}
            if(psize == 2)
                param = [param,0,1];
            elseif(psize == 3)
                param = [param,1];
            elseif(psize > 4 | psize == 1)
                error('Wrong number of parameter')
            end
            
            obj.mdl = mdl;
            obj.a = param(1);
            obj.range = param(2);
            obj.nugget = param(3);
            obj.smooth = param(4);
            obj.opt = opt;
            obj.link = link;
            
        case {'student'}
            if(psize == 2)
                param = [param,0,1];
            elseif(psize == 3)
                param = [param,1];
            elseif(psize > 4 | psize == 1)
                error('Wrong number of parameter')
            end
            
            obj.mdl = mdl;
            obj.df = param(1);
            obj.range = param(2);
            obj.nugget = param(3);
            obj.smooth = param(4);
            obj.opt = opt;
            obj.link = link;
                
        case {'t2'}
            if(psize == 3)
                param = [param,0,1];
            elseif(psize == 4)
                param = [param,1];
            elseif(psize > 5 | psize < 3)
                error('Wrong number of parameter')
            end
            
            obj.mdl = mdl;
            obj.df = param(1);
            obj.a = param(2);
            obj.range = param(3);
            obj.nugget = param(4);
            obj.smooth = param(5);
            obj.opt = opt;
            obj.link = link;
    end        
            
end %function

end % methods
end %classdef














