function [NI_META, AIN, DCode] = load_NI_data(NIFileName)
    
    NI_META=load_meta(sprintf('%s.meta', NIFileName));
    nFileBytes = NI_META.fileSizeBytes;
    nChan = NI_META.nSavedChans;
    nFileSamp = nFileBytes / (2 * nChan);
    fprintf('Load NI DATA\nn_channels: %d, n_file_samples: %d\n', nChan, nFileSamp);
    fprintf('Recording Last %04d seconds %03d mins\n', floor(nFileSamp./NI_META.niSampRate), floor(nFileSamp./NI_META.niSampRate/60));
    m = memmapfile(sprintf('%s.bin',NIFileName), 'Format', {'int16', [nChan, nFileSamp], 'x'},'Writable', false);
    NI_rawData = m.Data.x;

    fI2V = NI_META.niAiRangeMax/32768;
    MN = NI_META.snsMnMaXaDw(1);
    MA = NI_META.snsMnMaXaDw(2);
    XA = NI_META.snsMnMaXaDw(3);
    DW = NI_META.snsMnMaXaDw(4);
    
    digCh = MN + MA + XA + 1;
    
    AIN=double(NI_rawData(1,:))*fI2V;
    
%     digital0=NI_rawData(digCh,:);
% 
%     
%     CodeAll = diff(digital0);
%     DCode.CodeLoc = find(CodeAll>0);
%     DCode.CodeVal = CodeAll(DCode.CodeLoc);
%     

    digital0=NI_rawData(digCh,:);
    CodeAll = diff(digital0);
    DCode.CodeLoc = [find(CodeAll~=0)];
    digital0(1)=[];
    DCode.CodeVal = digital0(DCode.CodeLoc);
    DCode.CodeLoc = [nan, DCode.CodeLoc];
    DCode.CodeVal = [0,DCode.CodeVal];
    % this means SYNC code occurs with Onset Code, just fix it as simple as
    % possible, you can split 65 as 64+1 ans 63 as 64-1, but I will just ignore it now
%     DCode.CodeVal(DCode.CodeVal==63)=64;
%     DCode.CodeVal(DCode.CodeVal==65)=64;
%     DCode.CodeVal(DCode.CodeVal==3)=2;
    fprintf('Load Event Data\n')
    all_code = unique(DCode.CodeVal);
    for code_now = all_code
        fprintf('Event %d %d times\n', code_now, sum(DCode.CodeVal==code_now))
    end

    
    % Convert Data Into MS
    DCode.CodeTime = 1000*DCode.CodeLoc/NI_META.niSampRate;
    % ReSample
    [p, q] = rat(1000 / NI_META.niSampRate);
    AIN = resample(AIN, p, q);

end

