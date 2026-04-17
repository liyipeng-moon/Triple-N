function SI = CalcSI(x1,x2)
% selectivity (d-prime) between 2 arrays.
SI = (mean(x1) - mean(x2)) / sqrt((var(x1) + var(x2)) / 2);
