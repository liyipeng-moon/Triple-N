function bin_scatter(x,y,binx,biny)

density_plot = zeros([length(binx)-1, length(biny)-1]);
for xx = 1:length(binx)-1
    for yy = 1:length(biny)-1
        units_within = find(x>binx(xx) & x<=binx(xx+1) & y>biny(yy) & y<=biny(yy+1));
        density_plot(xx,yy) = length(units_within);
    end
end

density_plot = density_plot';
cm_here = colormap_matplotlib('coolwarm');
cm_here = cm_here(floor([257:513]/2),:);
imagesc(binx,biny,density_plot)
cm_here(1,:) = [1,1,1];


clim([0.9,max(density_plot(:))]);
set(gca,'ColorScale','log')
colorbar;
set(gca,'Colormap',cm_here)  
box on
set(gca,'YDir','normal')
plot([binx(1),binx(end)], [biny(1),biny(end)],'Color','k','LineStyle','--','LineWidth',1)
xlim([binx(1),binx(end)])
ylim([biny(1),biny(end)])
bin_stap = (binx(end)-binx(1))/5;
xticks([binx(1):bin_stap:binx(end)])
yticks([biny(1):bin_stap:biny(end)])
