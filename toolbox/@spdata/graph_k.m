function [xbreaks,xcdf,ybreaks,ycdf]=graph_k(obj,id,varargin)
%% Summary
%
% The function calculate the cdf of the probability integral 
% transformation
% 
%% CALL
%
%  [xbreaks,xcdf,ybreaks,ycdf]=graph_k(obj,id,varargin)
%
%  xbreaks: Point at which the empirical cdf is computed
%  xcdf   : Emprirical cdf
%  ybreals: Points at which the fitted cdf is computed
%  ycdf   : Fitted cdf approximated by a large number 
%           of simulation (default 2000). Use obj.spmodel.
%  id     : Identification number for the group of lag distance
%
% Additional options:
% 
%  (...,'plot')       : Overlay plot of xcdf et ycdf
%  (...,'plot_add')   : like 'plot' but no new figure is created
%  (...,'nsample', N) : Size N the sample used to approximate ycdf
%  (...,'ndraw', n)   : Number of points n to evaluate the fitted line
%  (...,'xch')        : Use the empirical cdf using exchengeability
%  (...,'cvm')        : Use the classical empirical cdf;
%  (...,'line_emp',S) : String S describing xcdf line (see function plot)
%  (...,'line_fit',S) : String S describing ycdf line (see function plot)
%  (...,'emp_only')   : Compute only xcdf


%%%%%%%%% HEADER %%%%%%%%%%%% 
do_plot = 0;
do_new = 0;
do_theo = 1;
nsample = 2000;
ndraw = 200;
line_emp = '-b';
line_fit = '-g';
method = 'xch';


for(ii = 1:numel(varargin))
    switch varargin{ii}
        case 'plot_add'
            do_plot = 1;
        case {'plot_new','plot'}
            do_plot = 1;
            do_new = 1;
        case 'nsample'
            nsample = varargin{ii+1};
        case 'ndraw'
            ndraw = varargin{ii+1};
        case 'line_emp'
            line_emp = varargin{ii+1};
        case 'line_fit'
            line_fit = varargin{ii+1};
        case {'normal','cvm'}
            method = 'cvm';
        case 'xch'
            method = 'xch';
        case 'emp_only'
            do_theo = 0;
    end   
end

%%%%%%%%% BODY %%%%%%%%%%%

x = copula.transform(obj.get_pairs(id), method);
xbreaks = unique(sort(x));
%xbreaks  = linspace(0,1)';
xcdf = ecdf_fast(x,xbreaks);
xcdf = xcdf(2:(end-1));

% theoritical curve approximate
if(do_theo)
    rho = obj.sp.rho(obj.mids);
    y = copula.transform(obj.sp.sim_biv(nsample,rho(id)));
    ybreaks = linspace(0,1,ndraw)';
    ycdf = ecdf_fast(y,ybreaks);
    ycdf = ycdf(2:(end-1));
end

if(do_plot)
    
    if(do_new) 
        figure; 
    else
        hold on
    end
    stairs (xbreaks,xcdf,line_emp)
    
    if(do_theo)
        hold on 
        plot(ybreaks,ycdf,line_fit)
    end
    
    if(do_new) hold off; end
end
