function start_parfor()
    if(isempty(gcp('nocreate')))
        numCores = feature('numcores');
        parpool(numCores-2);
    end
end