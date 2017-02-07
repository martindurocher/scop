function padj = padjust(pv,method)
%
%% DESCRIPTION
%
% Adjusting p-value for multiple testing
%
%% CALL 
%
%    padj = padjust(pv,method)
%
% pv     : P-values of the null hypotheses
% method : Method to calculate the correction (see NOTE)
%
%% NOTE
%
% The method are the following
% BH : Benjamini and Hochberg (1995) False discovery rate
% BY : Benjamini and Yekutieli (2001) False discovery rate - dependance
% HL : Holm (1979) familywise error - dependence
% HC : Hochberg (1988) 
%

%defaut values
if(nargin < 2)
    method = 'BH';
end

% pivot if vectors
[n,m] = size(pv);
if(n == 1 && m > 1 )
    pv = pv';
    pivot = 1;
    [n,m] = size(pv);
elseif(n == 1 && m == 1)
    padj = pv;
    return        
else
    pivot = 0;
end

% Determine the denominator of the correction
switch method
    case {'BH','fda'}
        a = n ./(1:n); 
    case {'HL','HC'}
        a = n - (1:n) + 1;
    case 'BY'
        l = sum(1 ./(1:n));
        a = n ./(1:n) * l;
end

padj = zeros(n,m);

for(j = 1:m)
    [pv_x,pv_i] = sort(pv(:,j));
    
    % adjusting the values
    ap = a' .* pv_x;
    
    % correcting the violation
    for(k = 1:(n-1))
        if(ap(k) > ap(k+1))
            switch method 
                case 'HL'
                    ap(k) = max(ap(1:k));
                otherwise
                   ap(k) = min(ap(k:n));
            end
        end
    end
    ap(n) = max(ap(n),ap(n-1));
    
    padj(pv_i,j) = ap;
end

padj = min(padj,1);

% if was a vector return a vector
if(pivot)
   padj = padj'; 
end
   



