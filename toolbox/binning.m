function lst = binning(x, breaks)
%% SUMMARY: 
%
%  Return a cell structure listing which indices fall into 
%      classes defined by breaks points
%
%% CALL 
%
%       lst = binning(x, breaks)
%
%  lst    : List of indices (Cell)
%  x      : Sample (Vector) 
%  breaks : List of distance defining (Vector)

p = max(size(breaks));
   
lst = cell(p-1,1);
for i=2:p
   lst{i-1} = find(x >= breaks(i-1) & x < breaks(i)); 
end
