cd C:\Users\moonl\Desktop\NNN\NNN_Data\Raw\H5FILES
all_saved_pro = dir('*info*');
tr = [];
deg =[];
IT_session = [1:70, 88];
for ii = IT_session
    load(all_saved_pro(ii).name)
    deg(ii) = [trial_ML(1).VariableChanges.img_degree_h];
    tr(ii)=length(meta_data.onset_time_ms);
end
mean(tr(IT_session))
std(tr(IT_session))

fprintf('ses num: %d\n', ii)
fprintf('trial num: %d(%d)\n', floor(mean(tr(IT_session))),floor(std(tr(IT_session))))

sum(deg==10)
sum(deg==12)
cd ..
%%
cd C:\Users\moonl\Desktop\NNN\NNN_Data\Raw\H5FILES
all_saved_pro = dir('*info*');
for ii = [1:90]
    load(all_saved_pro(ii).name)
    ff(ii) = [trial_ML(1).VariableChanges.fixation_window];
    ii
end
