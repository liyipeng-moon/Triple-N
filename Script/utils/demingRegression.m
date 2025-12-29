function [slope,CI, slopes]=demingRegression(x,y)
        n=length(x);
        lambda=1;
        alpha=0.01;
        sxx=sum(x.^2);
        syy=sum(y.^2);
        sxy=sum(x.*y);
        slope=(-lambda*sxx+syy+sqrt(4*lambda*sxy^2+(lambda*sxx-syy)^2))/(2*sxy);
        boot_times = 1000;
        slopes = zeros([1,boot_times]);
        for b = 1:boot_times
            idx = randsample(n, n, true);
            x_now = x(idx);
            y_now = y(idx);
            sxx=sum(x_now.^2);
            syy=sum(y_now.^2);
            sxy=sum(x_now.*y_now);
            slopes(b)=(-lambda*sxx+syy+sqrt(4*lambda*sxy^2+(lambda*sxx-syy)^2))/(2*sxy);
        end
        CI = quantile(slopes, [alpha/2, 1-alpha/2]);
end