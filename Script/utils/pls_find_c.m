function [r,r2,best_c] = pls_find_c(feature, response, cv_num)

response = zscore(response);
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(feature, response', 10,'cv',20);

[~,best_c]=min(MSE(2,2:end));
cv = cvpartition(size(response,2), 'KFold', cv_num);
pred_response = zeros(size(response));

for fold = 1:cv.NumTestSets
    trainIdx = training(cv, fold);
    testIdx = test(cv, fold);
    features_train = feature(trainIdx, :);
    rsp_train = response(trainIdx);
    features_test = feature(testIdx, :);
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(features_train, rsp_train', best_c);
    pred_response(testIdx) = [ones(size(features_test, 1), 1) features_test] * BETA;
end

r = corr(response', pred_response');
r2 = calculateR2(response',pred_response');

end