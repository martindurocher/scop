%% Getting started with the scop package

%% Introduction
%
% This package provides tools to carry out some aspect of analyzing 
% spatial data from a copula-based semiparametric framework. 
% In this approach the marginal are estimated by the empirical cdf 
% and the dependance is described by a copula model. The euclidien distance
% is the only distance implemented and the data are assumed isotropic. The 
% main interest of the package is to realize Goodness-of-fit test.

%% A simple exemple
%
% First let create object that represente a spatial copula model,
% i.e. a spmodel. Here a normal copula with 
% an exponential link function is create. 

spn = spmodel('normal','exp',[.20,.1],'r');

% the third terms is a vector of parameter, which will dependend of the
% copula. For the normal copula the parameter are related to the RANGE, the
% NUGGET effect and (if necessary) a third parameter is the SMOOTH parameter.
% The fourth argument is use to specify which parameter are variable. Like
% above, the RANGE = .15 and the NUGGET = .1, but as 'r' implies that only
% the range is a variable and the NUGGET is fixed.

% Similarely, here chi-square and Student copulas are defined with
% rational quadratic and spherical link function. Notice that the Student 
% copula here have option 'rn', which mean that both the range and the 
% Nugget effect are variables. 

spc = spmodel('chisq','rquad',[0,.4],'r');
spt = spmodel('student','sph', [5,.15],'r');

% Contrary to the normal copula, these two copula models have copula 
% specific parameters, that correspond to the first element of the parameter vector.
% vector. For instance the chi-square copula is centered (0) and the
% Student copula have 5 degree of freedom.

%%
% Let creates coordinates and simulates a datasets from these coordinates 
% using the spmodel of normal copula. The result is two simulated normal
% vector in N x 2 matrix, with coordinate chosen at random in the unit
% square. 

nsite = 300;
coord0 = rand(nsite,2);
data1 = spn.sim_coord(2, coord0)';

%% 
% Those data can be analyzed by using a spdata object, which joint the
% spatial data, the coordinate and a spmodel.

data0 = data1(:,1);
sdn = spdata(data0,coord0,spn);
sdc = spdata(data0,coord0,spc);
sdt = spdata(data0,coord0,spt);

% when creating the spdata object pseudo-observation in [0,1] are
% computed using the ranks. 
% If the margins are already in the [0,1] space,
% the 'direct' option can be specified to used them directly.

[data0(1:5),sdn.ranks(1:5)]

sdn = spdata(data0,coord0,spn,'direct');
[data0(1:5),sdn.ranks(1:5)]

%%
% The parameter of the model can be fitted as follow;

[sdn, boot_param] = sdn.fit('full','nboot',100);
sdn.sp %the updata spmodel object
prctile(boot_param,[5,50,95])'

sdc = sdc.fit('pairwise','htol',.1);
sdc.sp


%%
% In the case of the sdn object, the maximum likelihood is used 
% to estimate the parameter of the normal copula and return also 
% 100 bootstraps samples to obtained confident interval.
% In the case of the sdc object, the chi-square copula does not have an
% expression of the full likelihood and the pairwise likelihood is used
% instead. The parameter htol indicate a treshold where paired observation
% of distance further than htol are ignored.

%%
% Pairwise analysis of the bivariate copula arecan be performed among group 
% group of lag distance. The following function specify breaks points and 
% the utilization of the medium distance in each bin as middle points.
% Exploration graph are then provided to examine the evolution of the 
% measure of association and the kendall transformation in each bins. 

mybrks = [0,prctile(sdn.pdistance,(1:10)/2)]
sdn = sdn.set_breaks(mybrks,'med');

sdn.bin_size
sdn.plot_bin()
sdn.graph_k(1,'plot','line_fit','.r');

%%
% to perform a goodness of fit test, first set the group of
% lag distance. Then, perform bootstrap of the test statistics.

gof_out = sdn.test_gof('xchk',200); %very time consuming

%%
% The statistics and the pvalue can be obtained afterward.

[stat,pval] = gof_out.summary('xchk',0.05);
[stat;pval]'

% If there is K group of lag distance, the vector stat is of size 2K.
% The first K values are the test statistics for each groups and the
% following value are the cumulatives. Formal testing is obtained
% From the p-values (pval) of the cumulative.
%
% An alternative method consist to identify which group of 
% lag distance are not following the model by using 
% False discovery rate, which adjust the p-values for 
% multiple hypothesis testing

padjust(pval(1:5),'fda')


%%
% Here the parameter of the simulation are estimated using directly
% the spmodel object

%Calculating the ranks and the distances
data0 = rankit(data0);
dist0 = dist(coord0');

spn.fit('full',dist0,data0)

%pairwise likelihood
pdist0 = pdist(coord0);
pairs0 = pairwise(data0);
pid = (pdist0 < .05); 
spn.fit('pairwise',pdist0(pid),pairs0(pid,:))

% Direct computation pdf and cdf for the non centered chi-square copula 
% can be long, but can be speed-up by using a predefined table of key 
% components. Here an illustration using the spmodel.set_ha that construct
% the required table for 1000 of points

%without speed-up
spc.a = .5;
tic;
spc.fit('pairwise',pdist0(pid),pairs0(pid,:)) % few minutes
toc

%with speed-up
spc = spc.set_ha(1000);
tic;
spc.fit('pairwise',pdist0(pid),pairs0(pid,:))
toc

%%
% Further developement should bring kriging capability

pts = (1:30)/31;
grid = expand_grid(pts);
sim = sdn.sim_cond(50,grid);

pred = prctile(sim,[5,50,95]);
sd = std(sim)

contourf(pts,pts,reshape(pred(2,:),30,30))
hold on
scatter(coord0(:,1),coord0(:,2))

contourf(pts,pts,reshape(sd,30,30))
hold on
scatter(coord0(:,1),coord0(:,2))
hold off

