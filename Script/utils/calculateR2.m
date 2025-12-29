function R_squared = calculateR2(y_true, y_pred)
    y_true = y_true(:);
    y_pred = y_pred(:);
    SS_tot = sum((y_true - mean(y_true)).^2);
    SS_res = sum((y_true - y_pred).^2);
    R_squared = 1 - (SS_res / SS_tot);
end