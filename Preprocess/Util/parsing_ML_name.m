function [a,b] = parsing_ML_name(ml_name)
    split = find(ml_name=='_');
    a = ml_name(1:split(1)-1);
    b = ml_name(split(1)+1:split(2)-1);
end