function [META_DATA, DCode] = load_IMEC_data(NIFileName)
    META_DATA=load_meta(sprintf('%s.meta', NIFileName));
    nFileBytes = META_DATA.fileSizeBytes;
    nChan = META_DATA.nSavedChans;
    nFileSamp = nFileBytes / (2 * nChan);
    fprintf('Load IMEC DATA\nn_channels: %d, n_file_samples: %d\n', nChan, nFileSamp);
    fprintf('Recording Last %04d seconds %03d mins\n', floor(nFileSamp./META_DATA.imSampRate), floor(nFileSamp./META_DATA.imSampRate/60));
    m = memmapfile(sprintf('%s.bin',NIFileName), 'Format', {'int16', [nChan, nFileSamp], 'x'},'Writable', false);
    digital0 = m.Data.x(385,:);

    CodeAll = diff(digital0);
    DCode.CodeLoc = find(CodeAll>0);
    DCode.CodeVal = CodeAll(DCode.CodeLoc);
    fprintf('Load Event Data\n')
    all_code = unique(DCode.CodeVal);
    for code_now = all_code
        fprintf('Event %d %d times\n', code_now, sum(DCode.CodeVal==code_now))
    end

    % Convert Data Into MS
    DCode.CodeTime = 1000*DCode.CodeLoc/META_DATA.imSampRate;
end

