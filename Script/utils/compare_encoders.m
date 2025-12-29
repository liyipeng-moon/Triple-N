function [r,r2,best_numbers] = compare_encoders(xs,y,cv_num)
for ll = 1:length(xs)
    [r(ll),r2(ll),best_numbers(ll)] = pls_find_c(xs{ll}, y, cv_num);
end
end