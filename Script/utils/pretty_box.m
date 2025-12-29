function pretty_box(cm_here)

h = findobj(gca,'Tag','Box');
med = findobj(gca, 'Tag', 'Median');
hWhisker1 = findobj(gca, 'Tag', 'Upper Whisker');
hWhisker2 = findobj(gca, 'Tag', 'Lower Whisker');
for j = 1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),cm_here(j,:))
    scatter(med(j).XData(2),(med(j).YData(2)+med(j).YData(1))/2,12,'filled','MarkerEdgeColor','k',MarkerFaceColor='w')
    set(hWhisker1, 'LineStyle', '-');
    set(hWhisker1(j),'Color',[0.5,0.5,0.5])
    set(hWhisker1, 'LineWidth', 1);
    set(hWhisker2, 'LineStyle', '-');
    set(hWhisker2, 'LineWidth', 1);
    set(hWhisker2(j),'Color',[0.5,0.5,0.5])
end
caps = findobj(gca, 'Tag', 'Upper Adjacent Value');
delete(caps)
caps = findobj(gca, 'Tag', 'Lower Adjacent Value');
delete(caps)
