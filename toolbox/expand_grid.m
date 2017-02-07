function z = expand_grid(x,y)
if(nargin < 2); y=x;end
[a,b] = meshgrid(x,y);
z =[a(1:end);b(1:end)]';
