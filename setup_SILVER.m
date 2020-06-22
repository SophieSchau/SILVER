%% Download data
url = 'https://someurl.com/SILVER_data.tar.gz';
gunzip(url, 'SILVER_data');
untar('SILVER_data/SILVER_data.tar','SILVER_data');

%% Example 1 data
system('cp SILVER_data/example1/* examples/example1_ranges/.');

%% Example 2 data
system('cp SILVER_data/example2/* examples/example2_multiple_tempres/.');

%% Example 3 data
system('cp SILVER_data/example3/* examples/example3_ep_efficiency/.');

%% Example 4 data
system('cp SILVER_data/example4/* examples/example4_psf/.');

%% Example 5 data
system('cp -r SILVER_data/example5/* examples/example5_gfactor/.');

%% Example 6 data
system('cp SILVER_data/example6/* examples/example6_phantom_simu/.');