classdef spdata
%%
%% SUMMARY:
%
%    Perform spatial analysis of dependance structure with copula. 
% Primary focus on goodness-of-fit tests. 
%
%% CALL
%
%   obj = spdata(data, coord, sp, ...);

% data   : Sample of locations (Vector n x 1) 
% coord  : Euclidean Coordinates (Matrix n x 2)
% sp     : A spmodel, see spmodel.m

properties
    data % Matrix of data (nsite x 1)
    ranks % Ranks of the data matrix
    coord % Matrix of the coordinates (nsite x ncoord)
    pairs % Matrix of all pairs of ranks (npairs x 2)
    distance % Matrix of distance (nsite x nsite)
    pdistance % Distances in shape of the pairs (npairs x 1)
    nsite % number of sites

    sp % spmodel object

    breaks % Breaks that define bins boundaries 
    bins % Cell array containing the index all pairs for each bins
    nbin % Number of bins
    mids % Representative distance of the bins
    type % How the ranks are obtained either 'ranks' or 'direct'
    
    HISTSIZE = 500;
    GOF_N_APPROX = 500; %Number of simulation for approximation test_gof
    PROGRESS = 2;
        
    HTOL = 0;
end % end properties

methods

function obj = spdata(data,coord, sp, varargin)
    
    %%%%%%%%%%%%%%% header %%%%%%%%%%%%%%%%%%
    obj.type = 'ranks';
    
    for(ii = 1:numel(varargin))
    switch varargin{ii}
        case 'ranks'
            obj.type = 'ranks';
        case 'direct'
            obj.type = 'direct';
    end
    end
    
    %%%%%%%%%%% BODY %%%%%%%%%%%%%%%%%%%%
    
    obj.data = data;
    
    switch obj.type
        case 'ranks'
            obj.ranks = rankit(data);
        case 'direct'        
            if(min(data)<0 || max(data)>1)
                error('Some values are outside [0,1]')
            end
            obj.ranks = data;
    end
    
    obj.nsite = size(data,1);
    obj.pairs = pairwise(obj.ranks);
    obj.coord = coord;
    obj.distance = dist(coord');
    obj.pdistance = pdist(coord);
    obj.sp = sp;
    obj.HTOL = max(obj.pdistance);

end

end %method

end % class