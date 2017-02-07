function V = distance_ind(x,method,opt,param)
%
%% SUMMARY
%
% Calculate distance between an empirical copula and the product
% copula. 
%   
%% CALL
%
%  V = copula.distance(x,method,option,param)
%
% V      : Distance criterion (return)
% X      : sample 1 (observation in row)
% method : Method applied to obtained pseudo data (char)
%
%% Author Martin du Rocher <martin.durocher@uqtr.ca>

if(nargin <3); opt = 'd'; end;
if(nargin <4); param = 100; end;

% All values pass to C must be real number
if(~all_isfinite(x)); error('Non finite element in X'); end

%Must be 2D
if(size(x,2) ~= 2); error('X must be 2D'); end

% Transformation
switch opt
    case 'r'
        x = rankit(x);
    case {'hist','h'}
        method = [method,'h'];
end

switch method        
    case {'cvmk','ksk','cvmkt'}
        x = copula.transform(x);
        
    case {'xchk','xchkt'}
        x = copula.transform(x,'xch');

    otherwise
        x = sortrows(x);
end

% Calculation of the empirical copula and criterion
switch method
     case 'cvm'
        V = src.cvm2_ind_c(x(:,1),x(:,2));
                            
    case 'xch'
        V = src.xch2_ind_c(x(:,1),x(:,2));
        
    case {'cvmk','xchk'}
        x = max(1e-8,x - x .* log(x)); % avoid infinite
        V = src.cvm1_ind_c(x);
        
    case {'cvmkt','xchkt'}    
        V = src.cvm1_indt_c(x);
    
    case 'ks'
        breaks1 = unique(sort(x(:,1)))';
        breaks2 = unique(sort(x(:,2)))';
        cdf_x = src.emp_copula_c(x,breaks1,breaks2);
        cdf_ind = breaks1' * breaks2;
        V = max(abs(cdf_ind(1:end)-cdf_x(1:end)));
        
    case 'ksh'
        breaks1 =  (1:(param-1))/param;
        cdf_x = src.emp_copula_c(x,breaks1,breaks1);
        cdf_ind = breaks1' * breaks1;
        V = max(abs(cdf_ind(1:end)-cdf_x(1:end)));
        
    case 'ksk'
        breaks =  unique(sort(x));
        cdf_x = ecdf_fast(x,breaks);
        cdf_ind = breaks - breaks.*log(breaks);
        V = max(abs(cdf_x(2:end-1) - cdf_ind));
           
    otherwise 
       error('Wrong choice of method')
end

 