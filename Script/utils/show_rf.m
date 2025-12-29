function show_rf(fit_result, matrix)
    sf = 11/size(matrix,2);
    h = size(matrix,2)/2;
    [X, Y] = meshgrid(1:size(matrix,2), 1:size(matrix,1));
    X = (X-h-0.5)*sf;
    Y = (Y-h-0.5)*sf;
    
    x0 = fit_result.x0;
    y0 = fit_result.y0;
    sigma = fit_result.sigma;
    theta = linspace(0, 2*pi, 100);
    
    x_ellipse = x0 + sigma * cos(theta);
    y_ellipse = y0 + sigma * sin(theta);
    
    xx_array = linspace(-5.5,5.5,size(X,1));
    ll = [-5.5,5.5];

    imagesc(xx_array, xx_array, matrix);

    axis equal tight; xlim(ll); ylim(ll);
    plot(x_ellipse, y_ellipse, 'r', 'LineWidth', 1);  
    clim([prctile(matrix(:),5),prctile(matrix(:),100)])