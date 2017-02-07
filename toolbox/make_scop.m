function make_scop()
%% SUMMARY
% To automatize compilation of the mex file in the projet

cdir = pwd;
mydir = fileparts(which('make_scop'));

lst_file={'pairwise_c.c',...
    'kendall_transform_c.c', 'kendallx_transform_c.c',...
    'emp_copula_c.c', 'empx_copula_c.c',...
    'cvm1_c.c', 'cvm2_c.c', 'xch2_c.c','xch2_ind_c.c',...
    'cvm2_ind_c.c','cvm1_ind_c.c','cvm1_indt_c.c'};

clear mex;
if(ispc)
    cd([mydir,'\+src']);
else
    cd([mydir,'/+src']);
end

for(ii = 1:size(lst_file,2))
    ifile = lst_file{ii};
    disp(ifile)
    mex(ifile)
end

cd(cdir)

