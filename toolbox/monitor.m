function monitor(j,e,s)
%
%% DESCRIPTION
%
%  Print progress in the command line 
%
%% CALL
%    monitor(j,end,S)
%
%    j: Current step of the progess
%    end: Total number of steps in the progess
%    S: Upgrade progress every S steps
% 
%% Author: Martin du Rocher <martin.durocher@uqtr.ca>

    if(nargin <3); s = round(e/20); end
    
    s = round(s);
    if(j==e)
        fprintf('|\n');
    elseif(mod(j,s*5)==0)
       fprintf('|');
    elseif(mod(j,s)==0)
       fprintf('.');
    end
    
end