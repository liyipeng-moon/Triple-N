function [params, fit_result, gof] = gaussian2D_fit_symmetric(matrix)
sf = 11/size(matrix,1);
[X, Y] = meshgrid(1:size(matrix,2), 1:size(matrix,1)) ;
h = size(matrix,2)/2;
x_data = (X(:)-h-0.5) .* sf;
y_data = (Y(:)-h-0.5) .* sf;
z_data = matrix(:);
initial_guess = estimate_initial_parameters(matrix);
ft = fittype('A * exp(-((x-x0).^2 + (y-y0).^2)/(2*sigma^2)) + background',...
    'independent', {'x', 'y'},...
    'dependent', 'z',...
    'coefficients', {'A', 'x0', 'y0', 'sigma', 'background'});
opts = fitoptions(ft);
opts.StartPoint = initial_guess;
opts.Lower = [0, -7, -7, 0.0001, -max(abs(matrix(:))).^2];
opts.Upper = [max(abs(matrix(:))).^2, 7, 7, 6, max(abs(matrix(:))).^2];
[fit_result, gof] = fit([x_data, y_data], z_data, ft, opts);
params = [fit_result.A, fit_result.x0, fit_result.y0, fit_result.sigma, fit_result.background];
gof = gof.adjrsquare;
end
function initial_guess = estimate_initial_parameters(matrix)
[max_val, ~] = max(matrix(:));
background = mean(matrix(:));
sigma_estimate = 1;
[a,b] = find(matrix==max(matrix(:)));
a=a(1);b=b(1);
initial_guess = [max_val - background, a, b, sigma_estimate, background];
end