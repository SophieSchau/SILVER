%% Download data (basic)
system('tar -zxvf ./SILVER_data_basic.tar.gz');

%% Example 1 data
system('cp SILVER_data_basic/example1/* examples/example1_ranges/.');
system('rm SILVER_data_basic/example1/*');

%% Example 2 data
system('cp SILVER_data_basic/example2/* examples/example2_multiple_tempres/.');
system('rm SILVER_data_basic/example2/*');

%% Example 3 data
system('cp SILVER_data_basic/example3/* examples/example3_ep_efficiency/.');
system('rm SILVER_data_basic/example3/*');

%% Example 4 data
system('cp SILVER_data_basic/example4/* examples/example4_psf/.');
system('rm SILVER_data_basic/example4/*');

%% Example 5 data
system('cp -r SILVER_data_basic/example5/* examples/example5_gfactor/.');
system('rm -r SILVER_data_basic/example5/*');

%% Example 6 data
system('cp SILVER_data_basic/example6/* examples/example6_phantom_simu/.');
system('rm SILVER_data_basic/example6/*');

%% Precalculated SILVER optimizations
system('cp SILVER_data_basic/precalculated/* examples/precalculated/.');

%% Cleanup
system('rm -r SILVER_data_basic');


%% Example 7 (in-vivo) data
if exist('./SILVER_data_invivo.tar.gz', 'file')
    system('tar -zxvf ./SILVER_data_invivo.tar.gz');
    system('cp -r SILVER_data_invivo/example7/* examples/example7_invivo/.');
    system('rm -r SILVER_data_invivo');
end

