function tab = crosstab2(x,y)
%
%% DESCRIPTION
%
% Enhence crosstab for binary outcome. 
%
%% CALL
%
%  tab = crosstab2(x,y)
%
% x,y : Binary sample
%
%% Author: martin du rocher <martin.durocher@uqrt.ca>



% if a vector
[m,n] = size(x);
if(n==1 & m>1); x = x'; end
  
[m,n] = size(y);
if(n==1 & m>1); y = y'; end 

tab = zeros(3,3);

tab(1,1) = sum(x==0 & y==0);
tab(1,2) = sum(x==0 & y==1);
tab(2,1) = sum(x==1 & y==0);
tab(2,2) = sum(x==1 & y==1);
tab(3,:) = tab(1,:) + tab(2,:);
tab(:,3) = tab(:,1) + tab(:,2);