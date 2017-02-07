function out = sim_coord(obj,nobs,coord)
%% Method: sim_dist, sim_coord
%
% SUMMARY: Simulate a multivariate copula according to a matrix
%          of distance or coordinates.
%
% CALL: out = spmodel.sim_coord(nobs,coord)
%       out = spmodel.sim_dist(nobs,distance)
%
%   nobs : number of simulations (numeric)
%  coord : Coordinates of the sites (matrix n x 2)
%  dist  : matrix of distances (n x n)
%
% NOTE: see also sim_copula_multi.m and rsmith.m
%

switch obj.mdl
    case 'ind'
        out = rand(nobs,size(coord,1));
        
    case {'normal','chisq', 'student','t2'}   
        distance = dist(coord');
        sigma = obj.rho(distance); 
        tmp = obj.theta();
        out = copula.sim_multi(nobs,obj.mdl,sigma,tmp(1:2));
            
    case 'smith'
       
        if(strcmp(obj.opt,'iso'))
            sigma = [obj.cov11,0;0,obj.cov11];
        else
            [obj.cov11,obj.cov12;obj.cov12,obj.cov22];
        end 
        out = copula.sim_multi(nobs,'smith',coord,sigma);
end
    