%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The example data here were generated from ASL twix data from a 3T
%   Siemens Verio scanner acquired with a radial trajectory. To run example
%   7 on your own data you need to provide the raw kspace data in the same
%   format as generated here.
%
%   Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

S = [68, 153, 306];

file_Uniform = {...
    ['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-01-21/meas_MID195_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_Unif_68_FID47999.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-20/meas_MID221_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_Unif_68_FID50705.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-27/meas_MID243_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_Unif_68_FID50955.dat']; ... % 68 spokes/frame
    ['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-01-21/meas_MID196_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_Unif_153_FID48000.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-20/meas_MID223_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_Unif_153_FID50707.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-27/meas_MID244_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_Unif_153_FID50956.dat'];... % 153 spokes/frame
    ['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-01-21/meas_MID197_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_Unif_306_FID48001.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-20/meas_MID224_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_Unif_306_FID50708.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-27/meas_MID245_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_Unif_306_FID50957.dat']... % 306 spokes/frame
    };
file_GR = {...
    ['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-01-21/meas_MID198_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_GR_FID48002.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-20/meas_MID225_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_GR_FID50709.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-27/meas_MID242_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_GR_FID50954.dat'] ...
    };
file_SILVER = {...
    ['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-01-21/meas_MID199_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_SR_68_153_306_FID48003.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-20/meas_MID226_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_SR_68_153_306_FID50710.dat'],['/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Sedona_SILVER/2020-02-27/meas_MID240_GR_GRE_to_CAPIASL_CV_nce_angio_TD600_SR_68_153_306_FID50952.dat'] ...
    };

for subj = 1:length(file_SILVER)
    for n = 1:length(S)
        twix_obj = mapVBVD(file_Uniform{n,subj},'ignoreSeg');
        kdata_Uniform{n,subj} = reshape(permute(twix_obj.image{''},[1,3,5,4,2]),twix_obj.hdr.Config.RawCol, [], twix_obj.hdr.Config.NPhs,twix_obj.hdr.Config.NAveMeas,twix_obj.hdr.Config.NChaMeas);        
            
        spoke_order = [];
        NFrames = twix_obj.hdr.Config.NPhs;
        NRepeats = twix_obj.hdr.Config.NLin/twix_obj.hdr.Config.NSeg;
        NSpokes = twix_obj.hdr.Config.NLin;
        for repeat = 1:NRepeats
            spoke_order = cat(2, spoke_order, repeat:NRepeats:NSpokes);
        end
        [~, idx] = sort(spoke_order);
        kdata_Uniform{n,subj} = kdata_Uniform{n,subj}(:,idx,:,:,:);
        
        % Phase correction
        for ii = 1:size(kdata_Uniform{n,subj},2) % each line
            for jj = 1:size(kdata_Uniform{n,subj},3) % each frame
                for kk = 1:size(kdata_Uniform{n,subj},4) % each encoding
                    p=angle(mean(reshape(kdata_Uniform{n,subj}(:,ii,jj,kk,:).*conj(kdata_Uniform{n,subj}(:,ii,jj,1,:)),[],1)));
                    kdata_Uniform{n,subj}(:,ii,jj,kk,:)=kdata_Uniform{n,subj}(:,ii,jj,kk,:).*exp(-1j*p);
                end
            end
        end
        kdata_Uniform{n,subj} = reshape(kdata_Uniform{n,subj},[],twix_obj.hdr.Config.NLin*twix_obj.hdr.Config.NPhs/S(n),twix_obj.hdr.Config.NAveMeas,twix_obj.hdr.Config.NChaMeas); 
        

        twix_obj = mapVBVD(file_GR{subj},'ignoreSeg');
        kdata_GR{n,subj} = reshape(permute(twix_obj.image{''},[1,3,5,4,2]),twix_obj.hdr.Config.RawCol, [], twix_obj.hdr.Config.NPhs,twix_obj.hdr.Config.NAveMeas,twix_obj.hdr.Config.NChaMeas);        
        kdata_GR{n,subj} = kdata_GR{n,subj}(:,idx,:,:,:);
        % Phase correction
        for ii = 1:size(kdata_GR{n,subj},2) % each line
            for jj = 1:size(kdata_GR{n,subj},3) % each frame
                for kk = 1:size(kdata_GR{n,subj},4) % each encoding
                    p=angle(mean(reshape(kdata_GR{n,subj}(:,ii,jj,kk,:).*conj(kdata_GR{n,subj}(:,ii,jj,1,:)),[],1)));
                    kdata_GR{n,subj}(:,ii,jj,kk,:)=kdata_GR{n,subj}(:,ii,jj,kk,:).*exp(-1j*p);
                end
            end
        end
        kdata_GR{n,subj} = reshape(kdata_GR{n,subj},[],twix_obj.hdr.Config.NLin*twix_obj.hdr.Config.NPhs/S(n),twix_obj.hdr.Config.NAveMeas,twix_obj.hdr.Config.NChaMeas); 
        
        
        
        
        twix_obj = mapVBVD(file_SILVER{subj},'ignoreSeg');
        kdata_SILVER{n,subj} = reshape(permute(twix_obj.image{''},[1,3,5,4,2]),twix_obj.hdr.Config.RawCol, [], twix_obj.hdr.Config.NPhs,twix_obj.hdr.Config.NAveMeas,twix_obj.hdr.Config.NChaMeas);        
        kdata_SILVER{n,subj} = kdata_SILVER{n,subj}(:,idx,:,:,:);
        % Phase correction
        for ii = 1:size(kdata_SILVER{n,subj},2) % each line
            for jj = 1:size(kdata_SILVER{n,subj},3) % each frame
                for kk = 1:size(kdata_SILVER{n,subj},4) % each encoding
                    p=angle(mean(reshape(kdata_SILVER{n,subj}(:,ii,jj,kk,:).*conj(kdata_SILVER{n,subj}(:,ii,jj,1,:)),[],1)));
                    kdata_SILVER{n,subj}(:,ii,jj,kk,:)=kdata_SILVER{n,subj}(:,ii,jj,kk,:).*exp(-1j*p);
                end
            end
        end
        kdata_SILVER{n,subj} = reshape(kdata_SILVER{n,subj},[],twix_obj.hdr.Config.NLin*twix_obj.hdr.Config.NPhs/S(n),twix_obj.hdr.Config.NAveMeas,twix_obj.hdr.Config.NChaMeas); 
        
    end
end

SILVER_twix = twix_obj;

% Useful parameters
Mat_size = SILVER_twix.hdr.Config.BaseResolution;
NFrames = SILVER_twix.hdr.Meas.NPhs;
NRepeats = SILVER_twix.hdr.Config.NLin/SILVER_twix.hdr.Config.NSeg;
NSpokes = SILVER_twix.hdr.Config.NLin;
NSamps = SILVER_twix.hdr.Config.NColMeas;
NSubj = size(kdata_SILVER,2);



save('examples/example7_invivo/example_kdata_SILVER_68_153_306.mat', 'kdata_SILVER', '-v7.3');
save('examples/example7_invivo/example_kdata_GR_68_153_306.mat', 'kdata_GR', '-v7.3');
save('examples/example7_invivo/example_kdata_UNIFORM_68_153_306.mat', 'kdata_Uniform', '-v7.3');
save('examples/example7_invivo/example_params_68_153_306.mat', 'Mat_size', 'NFrames', 'NRepeats', 'NSpokes', 'NSamps', 'NSubj','S');