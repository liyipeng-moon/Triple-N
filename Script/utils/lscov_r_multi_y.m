function [r] = lscov_r_multi_y(feature, response, component_number)
response = zscore(response);
cv = cvpartition(size(response,1), 'KFold', 2);
trainIdx = training(cv, 1);
testIdx = test(cv, 1);

features_train = feature(trainIdx, :);
rsp_train = response(trainIdx,:);
rsp_test = response(testIdx,:);

pred_response = zeros(size(rsp_test));
for t = 1:size(response, 2)
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(features_train, rsp_train(:,t), component_number);
    features_test = feature(testIdx, :);
    pred_response(:,t) = [ones(size(features_test, 1), 1) features_test] * BETA;
end
r = corr(rsp_test, pred_response);
r = diag(r);
% figure;plot(diag(r).^2)
end
