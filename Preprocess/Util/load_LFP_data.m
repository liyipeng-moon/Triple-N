function [prb_meta, lfp_save] = load_LFP_data(local_path)
%     old ver
%     LF_file = dir([NIFileName,'\*lf.bin']);
%     LF_META_file = dir([NIFileName,'\*lf.meta']);
%     LFP_META=load_meta(fullfile(NIFileName,LF_META_file.name));
%     nFileBytes = LFP_META.fileSizeBytes;
%     nChan = LFP_META.nSavedChans;
%     nFileSamp = nFileBytes / (2 * nChan);
%     fprintf('Load LFP DATA\nn_channels: %d, n_file_samples: %d\n', nChan, nFileSamp);
%     fprintf('Recording Last %04d seconds %03d mins\n', floor(nFileSamp./LFP_META.imSampRate), floor(nFileSamp./LFP_META.imSampRate/60));
%     m = memmapfile(fullfile(NIFileName, LF_file.name), 'Format', {'int16', [nChan, nFileSamp], 'x'},'Writable', false);
%     NI_rawData = m.Data.x;
% 
%     fI2V = LFP_META.imAiRangeMax/32768;
%     LFP_IN=double(NI_rawData(1:end-1,:))*fI2V;
%     % Convert Data Into MS
%     % ReSample
%     
% 
%     example_channel = resample(LFP_IN(1,:), 1000, floor(LFP_META.imSampRate));
%     LFP_resampled = repmat(example_channel, [384,1]);
%     for cc = 1:size(LFP_IN,1)
%         LFP_resampled(cc,:) = resample(LFP_IN(cc,:), 1000, floor(LFP_META.imSampRate));
%     end

    cd(local_path)
    cd LFPprep\

    jsonText = fileread('binary.json');
    lfp_sr = jsondecode(jsonText).kwargs.sampling_frequency;
    
    jsonText = fileread('probe.json');
    prb_meta = jsondecode(jsonText).probes;
    
    filename = 'traces_cached_seg0.raw';
    mmap = memmapfile(filename,'Format','int16');
    data = mmap.Data;
    data = reshape(data, [384,length(data)./384]);
    data = double(data);
%     for c_idx = 1:384
%         data(c_idx,:) = zscore(data(c_idx,:));
%     end

    
    
    bad_channel = [];
    load ../channel_labels.mat
    for cc = 1:size(channel_labels,1)
        if(channel_labels(cc,1)~='g' || (mod(str2num(prb_meta.contact_ids{cc}(2:end)),384)==191))
            bad_channel = [bad_channel, cc];
        end
    end

    lfp_save = zeros([192, size(data,2)]);

    all_depth = prb_meta.contact_positions;
    depth_vals = unique(all_depth(:,2));
    for depth = 1:192
        this_channel = find(all_depth(:,2)==depth_vals(depth));
        if(~isempty(intersect(this_channel, bad_channel)))
            this_channel(this_channel==intersect(this_channel, bad_channel))=[];
        end
        lfp_save(depth, :) = mean(data(this_channel,:),1);
    end
    
    prb_meta.depth_vals = depth_vals;
    prb_meta.lfp_sr = lfp_sr;
    cd ..

    
end

