%% Download data
url_simu = './SILVER_data.tar.gz';
url_invivo = './SILVER_data_invivo.tar.gz';
gunzip(url_simu);
gunzip(url_invivo);
untar('SILVER_data.tar');
untar('SILVER_data_invivo.tar')

%% Example 1 data
system('cp SILVER_data/example1/* examples/example1_ranges/.');
system('rm SILVER_data/example1/*');

%% Example 2 data
system('cp SILVER_data/example2/* examples/example2_multiple_tempres/.');
system('rm SILVER_data/example2/*');

%% Example 3 data
system('cp SILVER_data/example3/* examples/example3_ep_efficiency/.');
system('rm SILVER_data/example3/*');

%% Example 4 data
system('cp SILVER_data/example4/* examples/example4_psf/.');
system('rm SILVER_data/example4/*');

%% Example 5 data
system('cp -r SILVER_data/example5/* examples/example5_gfactor/.');
system('rm -r SILVER_data/example5/*');

%% Example 6 data
system('cp SILVER_data/example6/* examples/example6_phantom_simu/.');
system('rm SILVER_data/example6/*');
% 
%% Example 7 data
system('cp -r SILVER_data_invivo/* examples/example7_invivo/.');
system('rm -r SILVER_data_invivo/*');

%% Precalculated SILVER optimizations
system('cp SILVER_data/precalculated/* examples/precalculated/.');

%% Final cleanup
system('rm SILVER_data*')
