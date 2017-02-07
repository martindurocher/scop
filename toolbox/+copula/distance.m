function V = distance(x,y,method,option,param)
%
%% SUMMARY
%
% Calculate distance between two empirical copula. 
%   
%% CALL
%
%  V = distance_copula(x,y,method,option,param)
%
% V      : Distance criterion (return)
% X      : sample 1 (observation in row)
% Y      : sample 2 (observation in row)
% method : Method applied to obtained pseudo data (char)
% option : option for the method (char)
% param  : optional parameter (numeric)
%
%% NOTE
%
% The available method are: Cramer Von Mises ()'cvm'), 
%  Kolmogorov-Smirnov 'ks' and Cramer Von Mises with 
%  exchangeability 'xch'
%
% There is 3 option. One is of using all sampled values ('emp') to 
%  calculate the empricical copula. The second is to use a predefined set of 
%  regular points ('hist'). In that case param is the number points.
%  The third option is to use a probability integral transformation ('kendall')
%  that create a univariate sample and use all sampled values
%
%% EXEMPLE 
% 
%  u1 = rand(100,2);
%  u2 = rand(150,2);
%
%  ans = copula.distance(u1,u2,'cvm','emp')
%  ans = copula.distance(u1,u2,'ks','hist',100)
%  ans = copula.distance(u1,u2,'xchk')
%  
%% Author Martin du Rocher <martin.durocher@uqtr.ca>

%Default parameters
if(nargin < 4); option = 'emp'; end;
if(nargin < 5); param = 100; end;

% All values pass to C must be real number
if(~all_isfinite(x)); error('Non finite element in X'); end
if(~all_isfinite(y)); error('Non finite element in Y'); end

%Must be 2D
if(size(x,2) ~= 2); error('X must be 2D'); end
if(size(y,2) ~= 2); error('Y must be 2D'); end

% for compatibility
switch method
    case 'cvmk'
        method = 'cvm';
        option = 'kendall';
    case 'xchk'
        method = 'xch';
        option = 'kendall';
end

% Transformation and ordering of the data
switch option
    case {'emp','hist'}
        x = sortrows(x);
        y = sortrows(y);
        
    case {'remp','rhist'}
        x = rankit(x);
        y = rankit(y);
        
        x = sortrows(x);
        y = sortrows(y);
        
    case 'kendall'
        if(strcmp(method,'xch'))
            x = copula.transform(x,'xch');
            y = copula.transform(y,'xch');
        else
            x = copula.transform(x);
            y = copula.transform(y);
        end
        
    otherwise
        error('Wrong option')
end

% Determine the points to evaluate the ecdf if necessary
switch option
    % Using a predefine set of points
    case {'rhist','hist'}
        breaks1    =  (1:(param-1))/param;
        breaks2    =  breaks1;
        
    %Using the union in case of kolmogorov-smirnov
    case {'remp','emp'}
        if(strcmp(method,'ks'))
            breaks1    =  unique(sort([x(:,1);y(:,1)]))';
            breaks2    =  unique(sort([x(:,2);y(:,2)]))';
        end
        
    case 'kendall'
        if(strcmp(method,'ks'))
            breaks =  unique(sort([x;y]));
        end
end

% Calculation of the empirical copula and distance criterion
method = [method,'_',option];
switch method
    
    case {'cvm_emp','cvm_remp'}
        V = src.cvm2_c(x(:,1),x(:,2),y(:,1),y(:,2));
        
    case {'cvm_hist','cvm_rhist'}
        cdf_x = src.emp_copula_c(x,breaks1,breaks2);
        cdf_y = src.emp_copula_c(y,breaks1,breaks2);
        V = mean(mean((cdf_x-cdf_y).^2));
        
    case 'cvm_kendall'
        V = src.cvm1_c(x(:,1),y(:,1));
        
    case {'ks_emp','ks_remp','ks_hist','ks_rhist'}           
        cdf_x = src.emp_copula_c(x,breaks1,breaks2);
        cdf_y = src.emp_copula_c(y,breaks1,breaks2);
        
        V = max(max(abs(cdf_x-cdf_y)));
        
    case 'ks_kendall'     
        cdf_x = ecdf_fast(x,breaks);
        cdf_y = ecdf_fast(y,breaks);
        V = max(abs(cdf_x-cdf_y));
        
    case {'xch_emp','xch_remp'}
        V = src.xch2_c(x(:,1),x(:,2),y(:,1),y(:,2));
        
    case {'xch_hist','xch_rhist'}
        cdf_x = src.empx_copula_c(x,breaks1,breaks2);
        cdf_y = src.empx_copula_c(y,breaks1,breaks2);
        V = mean(mean((cdf_x-cdf_y).^2));
        
    case 'xch_kendall'     
        V = src.cvm1_c(x(:,1),y(:,1));
                   
    otherwise 
       error('Wrong choice of method')
end


 