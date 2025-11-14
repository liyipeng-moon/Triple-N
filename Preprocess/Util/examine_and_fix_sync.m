function SyncLine = examine_and_fix_sync(DCode_NI, DCode_IMEC)


se = find(diff(bitand(DCode_NI.CodeVal,1))>0);
SyncLine.NI_time =  DCode_NI.CodeTime(1+se);
SyncLine.imec_time =  DCode_IMEC.CodeTime(DCode_IMEC.CodeVal==64);

fprintf('Syncing...\n')
fprintf('NI has %d edges while IMEC has %d\n', length(SyncLine.NI_time), length(SyncLine.imec_time))

d1 = diff(SyncLine.imec_time);
d1 = d1(2:end);
d2 = diff(SyncLine.NI_time);
d2 = d2(2:end);
figure;set(gcf,'Position',[67 100 325 900])
nexttile
plot(d1);ylim([950,2000])
title(sprintf('IMEC max diff is %f', max(d1)))
nexttile
plot(d2);ylim([950,2000])
title(sprintf('NI max diff is %f', max(d2)))
xlabel('# of rising edge')
ylabel('time lag between edges')


if(length(SyncLine.NI_time)~=length(SyncLine.imec_time))
    warning('Sync Fail! Fixing...\n')
    keyboard
    index = find(d2 > 1200);
    
    if ~isempty(index)
        for i = 1:length(index)
            idx = index(i)+i;
            new_value = (SyncLine.NI_time(idx) + SyncLine.NI_time(idx + 1)) / 2;
            SyncLine.NI_time = [SyncLine.NI_time(1:idx), new_value, SyncLine.NI_time(idx+1:end)];
        end
    end
    d1 = diff(SyncLine.imec_time);
    d1 = d1(2:end);
    d2 = diff(SyncLine.NI_time);
    d2 = d2(2:end);
    nexttile
    plot(d1);ylim([950,2000])
    title(sprintf('IMEC max diff is %f', max(d1)))
    nexttile
    plot(d2);ylim([950,2000])
    title(sprintf('NI max diff is %f', max(d2)))
    xlabel('# of rising edge')
    ylabel('time lag between edges')

else
    fprintf('Sync Success!\n')
end

terr = zeros([1,length(SyncLine.NI_time)]);
for ii = 1:length(SyncLine.NI_time)
    terr(ii)=SyncLine.NI_time(ii)-SyncLine.imec_time(ii);
end
nexttile
plot(terr)
ylim([-10,10])
end