function out = sim_dist(obj,nobs,distance)
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
        out = rand(nobs,size(distance,2));
        
    case {'normal','chisq','student','t2'}
        sample_sigma = obj.rho(distance);
        tmp = obj.theta();
        out = copula.sim_multi(nobs,obj.mdl,...
            sample_sigma,tmp(1:2));
    case 'smith'
        error('Not available')
        
end
    

