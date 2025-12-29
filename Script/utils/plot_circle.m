function plot_circle(fit_result,high_light)


    x0 = fit_result.x0;
    y0 = fit_result.y0;
    sigma = fit_result.sigma;
    theta = linspace(0, 2*pi, 100);
    x_ellipse = x0 + sigma * cos(theta);
    y_ellipse = y0 + sigma * sin(theta);

    scatter(x0,y0, 3,'filled','MarkerFaceColor',[0,0,0],MarkerFaceAlpha=0.5)
    plot(x_ellipse, y_ellipse, Color=[0 0 0 0.1], LineWidth=0.1);

    if(high_light)
        scatter(x0,y0,8,'filled','MarkerFaceColor',[1,0,0])
        plot(x_ellipse, y_ellipse, Color=[1 0 0 1], LineWidth=0.5);
    end