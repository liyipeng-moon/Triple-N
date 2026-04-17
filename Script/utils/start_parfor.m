function start_parfor()
% Initializes a parallel pool using most available CPU cores.
if(isempty(gcp('nocreate')))
    numCores = feature('numcores');
    parpool(numCores-2);
end
end