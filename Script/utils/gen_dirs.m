function [proc_dir,raw_dir] = gen_dirs(root_dir)
proc_dir = fullfile(root_dir,"NNN_Data/Processed/");
raw_dir = fullfile(root_dir,'NNN_Data/Raw/H5FILES/');