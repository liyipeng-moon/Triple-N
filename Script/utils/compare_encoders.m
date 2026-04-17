function [r,r2,best_numbers] = compare_encoders(xs,y,cv_num)
% Iterates over a set of input features xs to evaluate their predictive
% performance on y.
% 
for ll = 1:length(xs)
    [r(ll),r2(ll),best_numbers(ll)] = pls_find_c(xs{ll}, y, cv_num);
end
end