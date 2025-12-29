function bin_scatter_eye(x,y,binx,biny)
density_plot = zeros([length(binx)-1, length(biny)-1]);
for xx = 1:length(binx)-1
    for yy = 1:length(biny)-1
        units_within = find(x>binx(xx) & x<=binx(xx+1) & y>biny(yy) & y<=biny(yy+1));
        density_plot(xx,yy) = length(units_within);
    end
end
density_plot = density_plot';
cm_here = colormap_matplotlib('blues');
imagesc(binx,biny,density_plot)
clim([0,max(density_plot(:))]);
cc = colorbar;
cc.Label.String='Trial number';
set(gca,'Colormap',cm_here)  
box on
set(gca,'YDir','normal')
xlim([binx(1),binx(end)])
ylim([biny(1),biny(end)])
bin_stap = (binx(end)-binx(1))/5;
xticks([binx(1):bin_stap:binx(end)])
yticks([biny(1):bin_stap:biny(end)])

mu = [mean(x), mean(y)];
Sigma = cov(x,y);
theta = linspace(0,2*pi,200);
circle = [cos(theta); sin(theta)];
ellip = 2 * (chol(Sigma))' * circle;
ellip(1,:) = ellip(1,:) + mu(1);
ellip(2,:) = ellip(2,:) + mu(2);
plot(ellip(1,:),ellip(2,:),'k-','LineWidth',0.25)