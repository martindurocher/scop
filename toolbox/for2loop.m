function [jj,kk]=for2loop(ii,n)
%
%% DESCRIPTION
%
%  Help to transform a double for loop in one. Usefol for parfor context
%
%% CALL
%
%  Exemple of translation
%
%  niter = n*m;                    |  for(j=1:n)
%  parfor (ii =1:niter)            |     for(k=1:m)
%      [j,k] = for2loop(ii,n)      |        (some statements)  
%      (some statements)           |          ...
%       ...                        |     end    
%  end                             |  end
%                                  |  

%% Author: Martin du Rocher <martin.durocher@uqtr.ca>

jj = floor((ii-1)/n)+1;
kk = mod((ii-1),n)+1;