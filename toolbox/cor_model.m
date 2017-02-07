function y = cor_model(x, link, param)
%%
%% SUMMARY:
%
%   Function that returns correlation coefficient value conditional 
%       to covariable x
%% CALL
%
%    y = cor_model(x, link, param)
%
%  x     : Distance
%  model : Link function name (String)
%  param : Copula specific parameter (vector = [range,nugget,smooth])
%
%% NOTE
%    
%   List of correlation models
%   1. exponential ('exp')
%   2. gaussian ('gauss')
%   3. power exponential ('pow')
%   4. rational quadratic ('rquad')
%   5. Linear ('lin')
%   6. Spherical ('sph')
%   7. Circular ('circ')
%
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>

%case of an independent model
if(param(1) == 0)
    y = zeros(size(x));
    y(x==0) = 1;
    return
end

    % Computing the specified correlation model
    switch link
        case 'exp'
            y = exp(-3 * (abs(x) / param(1)));
        case 'gauss'
            y = exp(-3 * (abs(x) / param(1)) .^2);
        case 'pow'
            y = exp(-3 * (abs(x) / param(1)) .^param(3));
        case 'rquad'
            y = 1 ./ (1 + 19.*(x / param(1)) .^2 );
        case 'lin'
            y = 1-max(0,min(x/param(1),1));
        case 'sph'
            xa = x/param(1);
            xa(xa>1) = 1;
            y = 1 - 1.5*xa+.5*xa.^3;
        case 'circ'
            xa = x/param(1);
            xa(xa>1) = 1;
            y = 2/pi*(acos(xa)-xa.*sqrt(1-xa.*xa));
        otherwise
            error('Invalid correlation model')
    end
    
    y(y<1) = (1-param(2)) * y(y<1);

end


