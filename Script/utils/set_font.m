set(gcf,'color','w')
hText = findall(gcf, 'Type', 'text');
for i = 1:length(hText)
    set(hText(i), 'FontName', 'Arial');
    set(hText(i), 'FontSize', 9);
end
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 9;
box off