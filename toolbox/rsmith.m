function out = rsmith(nobs, coord, sigma)
% This function generates random fields for the 2d smith model

%   coord: the coordinates of the locations
%  center: the center of the compact set - here I use a square
%    edge: the length of the edge of the square
%    nObs: the number of observations to be generated
%    grid: Does coord specifies a grid?
%  nSites: the number of locations
%   covXX: the parameters of the bivariate normal density
%     ans: the generated random field */

%%
% Adapted by <martin.durocher@uqtr.ca>
%
% source:
%  @Manual{,
%    title = {SpatialExtremes: Modelling Spatial Extremes},
%    author = {Mathieu Ribatet},
%    year = {2015},
%    note = {R package version 2.0-2},
%    url = {https://CRAN.R-project.org/package=SpatialExtremes},
%  }

%%%%%%%%%%%%%%%%%%%
% HEAD
%%%%%%%%%%%%%%%%%%
nsites = size(coord,1);

%We first center the coordinates to avoid repetition of
%  unnecessary operations in the while loop */

coord_min = min(coord);
coord_max = max(coord);
edge = max(coord_max-coord_min); 
center = (coord_min + coord_max)/2;

for jj = 1:2
    coord(:,jj) = coord(:,jj)- center(jj);
end

% The compact set need to be inflated */

edge = edge + 6.92 * sqrt(max([sigma(1,1), sigma(2,2)]));
lebesgue = edge * edge;

%%%%%%%%%%%%%%%%%%%%%%
% BODY
%%%%%%%%%%%%%%%%%%%%%%
M_2_PI = 2/pi;
det = sigma(1,1) * sigma(2,2) - sigma(1,2) * sigma(1,2);
uBound = 1 / (M_2_PI * sqrt(det)); 
itwiceDet = 1 / (2 * det);
  
if ((det <= 0) || (sigma(1,1) <= 0))
    error('The covariance matrix is not semi-definite positive');
end
if(sigma(1,2) ~= sigma(2,1))
    error('The covariance matrix is not symmetrical');
end

% Simulation according to the Schlather methodology.
out = zeros(nobs,nsites);
kk=zeros(1,nobs);
for (ii = 1:nobs)
    poisson = 0;
    nKO = nsites;
    
    while(nKO)
	    % The stopping rule is reached when nKO = 0 i.e. when each site
	    %   satisfies the condition in Eq. (8) of Schlather (2002) */
        
        poisson = poisson + exprnd(10);
        ipoisson = 1 / poisson;
        thresh = uBound * ipoisson;

        %We simulate points uniformly in [-r/2, r/2]^2
        u1 = edge * (rand(1)-0.5);
        u2 = edge * (rand(1)-0.5);
      
        %nKO = nsites;
        %This is the bivariate normal density with 0 mean and
        % cov. matrix [sigma(1,1), sigma(1,2); sigma(1,2), sigma(2,2)]
       
        c1 = coord(:,1) - u1;
        c2 = coord(:,2) - u2;
        y = exp((-sigma(2,2) * c1.^2 + ...
                2 * sigma(1,2) * c1 .* c2 - ...
                sigma(1,1) * c2.^2) ...
                 * itwiceDet) * thresh;
	  
        out(ii,:) = max([y'; out(ii,:)]);
	  
        nKO = nsites - sum(thresh <= out(ii,:)); 
    end
    
end
  
% Lastly, we multiply by the Lebesgue measure of the dilated
%  compact set
out = out * lebesgue;

