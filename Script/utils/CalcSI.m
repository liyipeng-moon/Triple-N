function SI = CalcSI(x1,x2)
    SI = (mean(x1) - mean(x2)) / sqrt((var(x1) + var(x2)) / 2);
